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
DEFAULT_MIN_ABUNDANCE = 1.0e-1
DEFAULT_WAVELENGTH_MIN_ANGSTROM = 1.0
DEFAULT_WAVELENGTH_MAX_ANGSTROM = 250.0
DEFAULT_WAVELENGTH_STEP_ANGSTROM = 0.05


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


def compute_spectral_contribution(
    *,
    wavelength_min_angstrom: float,
    wavelength_max_angstrom: float,
    wavelength_step_angstrom: float,
    min_abundance: float,
    abundance: str,
    output_dir: Path,
) -> Path:
    wavelengths = fixed_step_grid(
        wavelength_min_angstrom,
        wavelength_max_angstrom,
        wavelength_step_angstrom,
    )
    densities = np.array([DEFAULT_DENSITY_CM3])
    temperatures = np.geomspace(1.0e4, 1.0e8, 201)
    component_values = np.full(
        (densities.size, temperatures.size, wavelengths.size, 4),
        np.nan,
        dtype=float,
    )

    print(
        f"Computing min_abundance={min_abundance:.1e} on "
        f"{wavelengths.size} wavelengths from {wavelengths[0]:g} to "
        f"{wavelengths[-1]:g} Angstrom with delta_lambda={wavelength_step_angstrom:g} Angstrom."
    )
    spectrum = None
    for density_index, density in enumerate(densities):
        spectrum = ch.spectrum(
            temperatures,
            density,
            wavelengths,
            filter=(chfilters.gaussian, DEFAULT_FILTER_WIDTH_ANGSTROM),
            em=None,
            doContinuum=True,
            minAbund=float(min_abundance),
            abundance=abundance,
            verbose=False,
        )
        component_values[density_index] = np.stack(
            [
                np.asarray(spectrum.FreeFree["intensity"], dtype=float),
                np.asarray(spectrum.FreeBound["intensity"], dtype=float),
                np.asarray(spectrum.LineSpectrum["intensity"], dtype=float),
                np.asarray(spectrum.TwoPhoton["intensity"], dtype=float),
            ],
            axis=-1,
        )

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
    included_elements = (element_abundances >= float(min_abundance)).astype(np.int8)
    abundance_info = chio.abundanceRead(spectrum.AbundanceName)
    abundance_reference = "".join(abundance_info["abundanceRef"]).strip()

    output_dir.mkdir(parents=True, exist_ok=True)
    grid_tag = (
        f"wavelength={wavelength_min_angstrom:g}-{wavelength_max_angstrom:g}"
        f"-dlambda={wavelength_step_angstrom:g}"
    )
    output_path = output_dir / (
        f"spectral-contribution.{grid_tag}.AbundanceName={spectrum.AbundanceName}"
        f"-min_abundance={min_abundance:.1e}.nc"
    )
    spectral_units = "erg cm3 s-1 sr-1 Angstrom-1"
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
            "min_abundance": float(min_abundance),
            "min_abundance_rule": "element_abundance >= min_abundance",
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
            "chiantipy_version": ChiantiPy.__version__,
            "chianti_database_version": chio.versionRead(),
        },
    )
    dataset.to_netcdf(output_path)
    print(f"Saved {output_path}")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute a CHIANTI spectral contribution function.")
    parser.add_argument("--wavelength-min", type=float, default=DEFAULT_WAVELENGTH_MIN_ANGSTROM)
    parser.add_argument("--wavelength-max", type=float, default=DEFAULT_WAVELENGTH_MAX_ANGSTROM)
    parser.add_argument("--wavelength-step", type=float, default=DEFAULT_WAVELENGTH_STEP_ANGSTROM)
    parser.add_argument("--min-abundance", type=float, default=DEFAULT_MIN_ABUNDANCE)
    parser.add_argument("--abundance", default=DEFAULT_ABUNDANCE)
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
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
