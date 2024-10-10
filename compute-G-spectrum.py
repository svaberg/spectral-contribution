import numpy as np
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters

abundance = 'test.abund'  # This is a local test abundance file
abundance = None          # This will make ChiPy use the default abundance file

wvls = np.linspace(1,180, 2001)             # Wavelengths in Angstroms
densities = np.array([1e10])                # Densities in cm^-3 (I think the units are cm^-3)
temperatures = np.geomspace(1e4, 1e8, 201)  # Temperatures in K

D, T, W = np.meshgrid(densities, temperatures, wvls, indexing = 'ij')
Z = np.zeros(D.shape + (5,))

np.save(f"grid-density.npy", D)
np.save(f"grid-temperature.npy", T)
np.save(f"grid-wavelength.npy", W)

min_abundances = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
for min_abundance in min_abundances:

    print(f"Now computing {min_abundance=:.1e}.")

    # Generating the intensities for the given temperatures, wavelengths and densities
    for i, density in enumerate(densities):
        spectrum = ch.spectrum(temperatures, 
                           density, 
                           wvls, 
                           filter=(chfilters.gaussian, 0.1), 
                           em = None, 
                           doContinuum=True, 
                           minAbund=float(min_abundance),
                           abundance=abundance,
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

    print(f"Using abundance file {spectrum.AbundanceName}.")
    save_name = f"G_lambda_T-spectrum.AbundanceName={spectrum.AbundanceName}-{min_abundance=:2.1e}.npy"
    print(f"Saving file named {save_name=}")
    np.save(save_name, Z)

    if np.min(spectrum.Abundance[spectrum.Abundance>0]) >= min_abundance:
        print('Reached minimum abundance')
        break



