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

    fig.colorbar(im,
                 ax=[axs['C'], axs['E']],
                 label = r'Intensity (units TBD)',
                 extend = 'both',
                 location='right',)

    axs['A'].set_xlabel(r'Wavelength $\lambda$ (Å)')
    axs['A'].set_ylabel(r'Electron temperature $T_\text{e}$ (K)')
    axs['A'].set_title('Spectral contribution function')
    save_name = Path(file_name).stem
    fig.suptitle(f"File: {file_name}")
    plt.savefig(f"plot-components-{save_name}.png")
    plt.show()
    plt.close()


if __name__ == '__main__':
    import sys
    file = sys.argv[1]
    main(file)
