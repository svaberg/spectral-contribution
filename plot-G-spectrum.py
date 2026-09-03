from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import LogNorm


def main(file_name):
    dataset = xr.load_dataset(file_name)

    total = dataset["freefree"] + dataset["freebound"] + dataset["line"] + dataset["twophoton"]
    components = [total, dataset["freefree"], dataset["freebound"], dataset["line"], dataset["twophoton"]]
    titles = ["total", "freefree", "freebound", "line", "twophoton"]

    fig, axs = plt.subplot_mosaic("""AABC
    AADE""", figsize=(12,6))

    d_id = 0
    wavelength = dataset["wavelength"].values
    temperature = dataset["temperature"].values

    cmap = matplotlib.colormaps["plasma"].copy()
    cmap.set_bad('grey')
    cmap.set_extremes(under='black', over='white')
    cmap.set_under(cmap(0))

    norm = LogNorm(vmax=1e-24, vmin=1e-29)

    for component, title, (k, ax) in zip(components, titles, axs.items()):
        im = ax.pcolormesh(
            wavelength,
            temperature,
            component.isel(density=d_id).values,
            norm=norm,
            cmap=cmap,
            rasterized=True,
        )

        ax.set_yscale('log')
        ax.set_title(f"{title} term")

    fig.colorbar(
        im,
        ax=[axs['C'], axs['E']],
        label=r'$G_\lambda(T_\mathrm{e})$ (erg cm$^3$ s$^{-1}$ sr$^{-1}$ $\AA^{-1}$)',
        extend='both',
        location='right',
    )

    axs['A'].set_xlabel(r'Wavelength $\lambda$ (Å)')
    axs['A'].set_ylabel(r'Electron temperature $T_\text{e}$ (K)')
    axs['A'].set_title('Spectral contribution function')
    source_path = Path(file_name)
    save_name = source_path.stem
    fig.suptitle(
        rf"$\lambda={wavelength[0]:g}$--${wavelength[-1]:g}\,\AA$, "
        rf"$\Delta\lambda={dataset.attrs['wavelength_step']:g}\,\AA$, "
        rf"minimum abundance $={dataset.attrs['min_abundance']:.1e}$"
    )
    output_path = source_path.parent / f"plot-components-{save_name}.png"
    plt.savefig(output_path)
    print(f"Saved {output_path}")
    plt.close()


if __name__ == '__main__':
    import sys
    file = sys.argv[1]
    main(file)
