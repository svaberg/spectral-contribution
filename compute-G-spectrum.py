import numpy as np
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters
import xarray as xr

abundance = 'test.abund'  # This is a local test abundance file
abundance = None          # This will make ChiPy use the default abundance file

wvls = np.linspace(1,180, 2001)             # Wavelengths in Angstroms
densities = np.array([1e10])                # Densities in cm^-3 (I think the units are cm^-3)
temperatures = np.geomspace(1e4, 1e8, 201)  # Temperatures in K

Z = np.zeros((len(densities), len(temperatures), len(wvls), 4))

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

        # Z persists across decreasing min_abundance runs. ChiantiPy can
        # return NaN values at lower minAbund settings, so we skip those
        # entries and keep the existing finite value from an earlier pass.
        for idx in range(intensity.shape[0]):
            for jdx in range(intensity.shape[1]):
                if np.isnan(intensity[idx, jdx]):
                    continue
                else:
                    Z[i, idx, jdx, 0] = freefree[idx, jdx]
                    Z[i, idx, jdx, 1] = freebound[idx, jdx]
                    Z[i, idx, jdx, 2] = line[idx, jdx]
                    Z[i, idx, jdx, 3] = twophoton[idx, jdx]

    print(f"Using abundance file {spectrum.AbundanceName}.")
    save_name = f"G_lambda_T-spectrum.AbundanceName={spectrum.AbundanceName}-{min_abundance=:2.1e}.nc"
    print(f"Saving file named {save_name=}")
    xr.Dataset(
        data_vars={
            "freefree": (("density", "temperature", "wavelength"), Z[:, :, :, 0]),
            "freebound": (("density", "temperature", "wavelength"), Z[:, :, :, 1]),
            "line": (("density", "temperature", "wavelength"), Z[:, :, :, 2]),
            "twophoton": (("density", "temperature", "wavelength"), Z[:, :, :, 3]),
        },
        coords={
            "density": ("density", densities, {"units": "cm^-3"}),
            "temperature": ("temperature", temperatures, {"units": "K"}),
            "wavelength": ("wavelength", wvls, {"units": "Angstrom"}),
        },
        attrs={
            "abundance_name": spectrum.AbundanceName,
            "min_abundance": float(min_abundance),
            "filter_width": 0.1,
            "filter_width_units": "Angstrom",
        },
    ).to_netcdf(save_name)

    if np.min(spectrum.Abundance[spectrum.Abundance>0]) >= min_abundance:
        print('Reached minimum abundance')
        break
