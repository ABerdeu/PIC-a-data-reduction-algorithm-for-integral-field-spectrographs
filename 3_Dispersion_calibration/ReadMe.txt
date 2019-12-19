3_Dispersion_calibration -> Dispersion model calibration by joint estimate on the wavelength calibration and spectral calibration files (Appendices B and C)

- executable -> run the first estimate of the dispersion calibration according to the parameters set in 'procedures/parameters.m' (mainly set the section '%% Data', the other sections are general parameters that are generally untouched)

- refine_pos_dif_laws_parallel -> run the refinement of the dispersion parameters obtained after the execution of 'executable.m' (mainly set the variable %% Parameters -> path_cal, the others are general parameters that are generally untouched)

- analyze.m -> plot the analysis of the different iteration of the refinement obtained after the executaion of 'refine_pos_dif_laws_parallel .m' (set the parameters in the section %% Set parameters)