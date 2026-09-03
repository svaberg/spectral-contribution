from __future__ import annotations

import argparse
from pathlib import Path

import ChiantiPy
import ChiantiPy.core as ch
import ChiantiPy.tools.constants as chconstants
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.io as chio
import numpy as np
import xarray as xr


DEFAULT_ABUNDANCE = "sun_photospheric_2015_scott"
DEFAULT_DENSITY_CM3 = 1.0e10
DEFAULT_FILTER_WIDTH_ANGSTROM = 0.1
DEFAULT_MIN_ABUNDANCE = 1.0e-7
DEFAULT_WAVELENGTH_MIN_ANGSTROM = 1.0
DEFAULT_WAVELENGTH_MAX_ANGSTROM = 250.0
DEFAULT_WAVELENGTH_STEP_ANGSTROM = 0.05
STANDARD_MIN_ABUNDANCE_RAMP = (
    1.0e0,
    1.0e-1,
    1.0e-2,
    1.0e-3,
    1.0e-4,
    1.0e-5,
    1.0e-6,
    1.0e-7,
)


def fixed_step_grid(start: float, stop: float, step: float) -> np.ndarray:
    if start <= 0.0:
        raise ValueError("The minimum wavelength must be positive")
    if stop <= start:
        raise ValueError("The maximum wavelength must exceed the minimum wavelength")
    if step <= 0.0:
        raise ValueError("The wavelength step must be positive")

    interval_count_float = (stop - start) / step
    interval_count = int(round(interval_count_float))
    if not np.isclose(interval_count_float, interval_count, rtol=0.0, atol=1.0e-10):
        raise ValueError(
            f"Wavelength interval {start:g}--{stop:g} Angstrom is not divisible "
            f"by the requested step {step:g} Angstrom"
        )
    return np.linspace(start, stop, interval_count + 1)


def min_abundance_ramp(final_threshold: float) -> tuple[float, ...]:
    if not 0.0 < final_threshold <= 1.0:
        raise ValueError("The final minimum-abundance threshold must be in (0, 1]")

    thresholds = [
        threshold
        for threshold in STANDARD_MIN_ABUNDANCE_RAMP
        if threshold >= final_threshold
    ]
    if not np.isclose(thresholds[-1], final_threshold, rtol=1.0e-12, atol=0.0):
        thresholds.append(float(final_threshold))
    return tuple(thresholds)


def compute_spectral_contribution(
    *,
    wavelength_min_angstrom: float,
    wavelength_max_angstrom: float,
    wavelength_step_angstrom: float,
    min_abundance: float,
    abundance: str,
    processes: int,
    output_dir: Path,
) -> list[Path]:
    if processes < 1:
        raise ValueError("The process count must be at least one")

    wavelengths = fixed_step_grid(
        wavelength_min_angstrom,
        wavelength_max_angstrom,
        wavelength_step_angstrom,
    )
    densities = np.array([DEFAULT_DENSITY_CM3])
    temperatures = np.geomspace(1.0e4, 1.0e8, 201)
    # This array deliberately persists across the decreasing-abundance ramp.
    # If a later ChiantiPy calculation returns NaN at a grid cell, retain the
    # finite component values calculated at the preceding threshold.
    component_values = np.zeros(
        (densities.size, temperatures.size, wavelengths.size, 4),
        dtype=float,
    )

    spectrum_class = ch.spectrum if processes == 1 else ch.mspectrum
    spectrum_kwargs = {} if processes == 1 else {"proc": processes}
    spectrum_backend = "spectrum" if processes == 1 else "mspectrum"
    thresholds = min_abundance_ramp(min_abundance)
    threshold_text = ", ".join(f"{threshold:.1e}" for threshold in thresholds)
    print(
        f"Computing minimum-abundance ramp [{threshold_text}] on "
        f"{wavelengths.size} wavelengths from {wavelengths[0]:g} to "
        f"{wavelengths[-1]:g} Angstrom with delta_lambda={wavelength_step_angstrom:g} Angstrom "
        f"using {spectrum_backend} with {processes} process(es)."
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    grid_tag = (
        f"wavelength={wavelength_min_angstrom:g}-{wavelength_max_angstrom:g}"
        f"-dlambda={wavelength_step_angstrom:g}"
    )
    spectral_units = "erg cm3 s-1 sr-1 Angstrom-1"
    output_paths = []

    for ramp_index, threshold in enumerate(thresholds):
        print(
            f"Ramp step {ramp_index + 1}/{len(thresholds)}: "
            f"min_abundance={threshold:.1e}."
        )
        spectrum = None
        nan_cell_count = 0
        for density_index, density in enumerate(densities):
            spectrum = spectrum_class(
                temperatures,
                density,
                wavelengths,
                filter=(chfilters.gaussian, DEFAULT_FILTER_WIDTH_ANGSTROM),
                em=None,
                doContinuum=True,
                minAbund=float(threshold),
                abundance=abundance,
                verbose=False,
                **spectrum_kwargs,
            )
            total_intensity = np.asarray(spectrum.Spectrum["intensity"], dtype=float)
            calculated_components = np.stack(
                [
                    np.asarray(spectrum.FreeFree["intensity"], dtype=float),
                    np.asarray(spectrum.FreeBound["intensity"], dtype=float),
                    np.asarray(spectrum.LineSpectrum["intensity"], dtype=float),
                    np.asarray(spectrum.TwoPhoton["intensity"], dtype=float),
                ],
                axis=-1,
            )
            valid_cells = ~np.isnan(total_intensity)
            nan_cell_count += int((~valid_cells).sum())
            component_values[density_index][valid_cells] = calculated_components[valid_cells]

        if spectrum is None:
            raise RuntimeError("No spectra were calculated")

        # ChiantiPy's abundance gate currently considers elements H--Zn (Z=1--30).
        # Embed the values that actually control that gate so the output remains
        # reproducible even if the named CHIANTI abundance file later changes.
        atomic_numbers = np.arange(1, 31, dtype=np.int16)
        element_symbols = np.asarray(
            [symbol.capitalize() for symbol in chconstants.El[: atomic_numbers.size]],
            dtype="U2",
        )
        element_abundances = np.asarray(
            spectrum.Abundance[: atomic_numbers.size],
            dtype=float,
        )
        included_elements = (element_abundances >= float(threshold)).astype(np.int8)
        abundance_info = chio.abundanceRead(spectrum.AbundanceName)
        abundance_reference = "".join(abundance_info["abundanceRef"]).strip()

        output_path = output_dir / (
            f"spectral-contribution.{grid_tag}.AbundanceName={spectrum.AbundanceName}"
            f"-min_abundance={threshold:.1e}.nc"
        )
        dataset = xr.Dataset(
            data_vars={
                "freefree": (
                    ("density", "temperature", "wavelength"),
                    component_values[:, :, :, 0],
                    {"units": spectral_units},
                ),
                "freebound": (
                    ("density", "temperature", "wavelength"),
                    component_values[:, :, :, 1],
                    {"units": spectral_units},
                ),
                "line": (
                    ("density", "temperature", "wavelength"),
                    component_values[:, :, :, 2],
                    {"units": spectral_units},
                ),
                "twophoton": (
                    ("density", "temperature", "wavelength"),
                    component_values[:, :, :, 3],
                    {"units": spectral_units},
                ),
                "element_symbol": (
                    "atomic_number",
                    element_symbols,
                    {"long_name": "chemical element symbol"},
                ),
                "element_abundance": (
                    "atomic_number",
                    element_abundances,
                    {
                        "long_name": "elemental number abundance relative to hydrogen",
                        "units": "1",
                    },
                ),
                "element_included": (
                    "atomic_number",
                    included_elements,
                    {
                        "long_name": "element passes the ChiantiPy minAbund gate",
                        "flag_values": np.asarray([0, 1], dtype=np.int8),
                        "flag_meanings": "excluded included",
                    },
                ),
            },
            coords={
                "density": ("density", densities, {"units": "cm-3"}),
                "temperature": ("temperature", temperatures, {"units": "K"}),
                "wavelength": ("wavelength", wavelengths, {"units": "Angstrom"}),
                "atomic_number": (
                    "atomic_number",
                    atomic_numbers,
                    {"long_name": "atomic number"},
                ),
            },
            attrs={
                "description": "CHIANTI spectral contribution function G_lambda(T)",
                "abundance_name": spectrum.AbundanceName,
                "abundance_reference": abundance_reference,
                "abundance_scale": "number abundance relative to hydrogen (H = 1)",
                "min_abundance": float(threshold),
                "min_abundance_rule": "element_abundance >= min_abundance",
                "ramp_final_min_abundance": float(min_abundance),
                "ramp_step_index": ramp_index,
                "ramp_step_count": len(thresholds),
                "ramp_thresholds": threshold_text,
                "nan_cells_preserved_from_previous_step": nan_cell_count,
                "nan_fallback_rule": (
                    "Where the current total spectrum is NaN, retain component values "
                    "from the preceding minimum-abundance step"
                ),
                "filter_name": "gaussian",
                "filter_width": DEFAULT_FILTER_WIDTH_ANGSTROM,
                "filter_width_units": "Angstrom",
                "filter_width_definition": "Gaussian standard deviation",
                "wavelength_min": float(wavelengths[0]),
                "wavelength_max": float(wavelengths[-1]),
                "wavelength_step": float(wavelength_step_angstrom),
                "wavelength_units": "Angstrom",
                "wavelength_grid_spacing": "linear",
                "temperature_grid_spacing": "geometric",
                "ionization_equilibrium_name": spectrum.Defaults["ioneqfile"],
                "flux_type": spectrum.Defaults["flux"],
                "do_continuum": 1,
                "do_lines": 1,
                "spectrum_backend": spectrum_backend,
                "process_count_requested": processes,
                "chiantipy_version": ChiantiPy.__version__,
                "chianti_database_version": chio.versionRead(),
            },
        )
        dataset.to_netcdf(output_path)
        print(
            f"Saved {output_path} "
            f"({nan_cell_count} NaN cell(s) retained from the preceding ramp state)."
        )
        output_paths.append(output_path)

    return output_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute a CHIANTI spectral contribution function.")
    parser.add_argument("--wavelength-min", type=float, default=DEFAULT_WAVELENGTH_MIN_ANGSTROM)
    parser.add_argument("--wavelength-max", type=float, default=DEFAULT_WAVELENGTH_MAX_ANGSTROM)
    parser.add_argument("--wavelength-step", type=float, default=DEFAULT_WAVELENGTH_STEP_ANGSTROM)
    parser.add_argument(
        "--min-abundance",
        type=float,
        default=DEFAULT_MIN_ABUNDANCE,
        help="Final (lowest) threshold in the decreasing minimum-abundance ramp.",
    )
    parser.add_argument("--abundance", default=DEFAULT_ABUNDANCE)
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Use spectrum with 1 process (default), or mspectrum with N > 1 processes.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compute_spectral_contribution(
        wavelength_min_angstrom=args.wavelength_min,
        wavelength_max_angstrom=args.wavelength_max,
        wavelength_step_angstrom=args.wavelength_step,
        min_abundance=args.min_abundance,
        abundance=args.abundance,
        processes=args.processes,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
