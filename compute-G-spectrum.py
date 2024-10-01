import numpy as np

import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters


wvls = np.linspace(1,180, 2001)
densities = np.array([1e10])
temperatures = np.geomspace(1e4, 1e8, 201)
D, T, W = np.meshgrid(densities, temperatures, wvls, indexing = 'ij')
Z = np.zeros(D.shape + (5,))

np.save(f"D.npy", D)
np.save(f"T.npy", T)
np.save(f"W.npy", W)

# Setting the abundance and resolution. Meshgrid fixes the dimension differences

min_abundances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

for min_abundance in min_abundances:

    save_name = f"{min_abundance:2.1e}-G"
    print(f"{save_name=}")
    print(f'{min_abundance=}')

    # Generating the intensities for the given temperatures, wavelengths and densities
    for i, density in enumerate(densities):
        spectrum = ch.spectrum(temperatures, 
                           density, 
                           wvls, 
                           filter=(chfilters.gaussian, 0.1), 
                           em = None, 
                           doContinuum=True, 
                           minAbund=float(min_abundance),
                           abundance='sun_coronal_2012_schmelz_ext.abund', 
                        #    proc=1,
                           verbose=False,
                           )
        
        intensity = spectrum.Spectrum['intensity']
        freefree = spectrum.FreeFree['intensity']
        freebound = spectrum.FreeBound['intensity']
        line = spectrum.LineSpectrum['intensity']
        twophoton = spectrum.TwoPhoton['intensity']

        for idx in range(intensity.shape[0]):
            for jdx in range(intensity.shape[1]):
                if np.isnan(intensity[idx, jdx]):
                    continue
                else:
                    Z[i, idx, jdx, 0] = intensity[idx, jdx]
                    Z[i, idx, jdx, 1] = freefree[idx, jdx]
                    Z[i, idx, jdx, 2] = freebound[idx, jdx]
                    Z[i, idx, jdx, 3] = line[idx, jdx]
                    Z[i, idx, jdx, 4] = twophoton[idx, jdx]

    np.save(f"Z-{save_name}.npy",Z)
