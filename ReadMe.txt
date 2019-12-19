- 1_Data -> Convert the file names of the raw data downloaded from the data center to file names that are interpretable by the different procedures

- 2_Preprocessing_flat_and_dark -> Sensor calibration: dark, sky background, flat, noise model

- 3_Dispersion_calibration -> Dispersion model calibration by joint estimate on the wavelength calibration and spectral calibration files

- 4_PIC_operators_declaration -> Build the operators of the forward model and performs their autocalibration on the science data.

- 5_Inverse_problem_reduction -> Data reduction by inverse problem approach with the PIC forward model

- Add-ons -> general toolboxes needed for the different procedures