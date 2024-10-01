import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

D = np.load(f"D.npy")  # Density grid
T = np.load(f"T.npy")  # Temperature grid
W = np.load(f"W.npy")  # Wavelength grid


def main(file_name):
    Z = np.load(file_name, allow_pickle=True)
    
    fig, axs = plt.subplot_mosaic("""AABC
    AADE""", figsize=(12,6))
    
    d_id = 0
    
    cmap = matplotlib.colormaps["plasma"].copy()
    cmap.set_bad('grey')
    cmap.set_extremes(under='black', over='white')
    cmap.set_under(cmap(0))
   
    norm = LogNorm(vmax=1e-24, vmin=1e-29)

    titles = "intensity freefree freebound line twophoton".split()
    
    for i, (k, ax) in enumerate(axs.items()):

        im = ax.pcolormesh(W[d_id, ...], 
                        T[d_id, ...], 
                        Z[d_id, ..., i],
                        norm=norm,
                        cmap=cmap, 
                        rasterized=True,
                        )
        
        ax.set_yscale('log')
        ax.set_title(f"{titles[i]} term")
    
    fig.colorbar(im, 
                 ax=[axs['C'], axs['E']],
                # cax=cax, 
                label = r'Intensity (units TBD)', 
                extend = 'both',
                location='right',)

    axs['A'].set_xlabel('Wavelength $\lambda$ (Å)')
    axs['A'].set_ylabel(r'Electron temperature $T_\text{e}$ (K)')
    axs['A'].set_title(r'Contribution function $G(T_\text{e}, \lambda)$')
    save_name = file_name.split('/')[-1]
    fig.suptitle(f"File: {file_name}")
    plt.savefig(f"plot-components-{save_name}.png")
    plt.show()
    plt.close()

if __name__ == '__main__':
    import sys
    file = sys.argv[1]
    main(file)