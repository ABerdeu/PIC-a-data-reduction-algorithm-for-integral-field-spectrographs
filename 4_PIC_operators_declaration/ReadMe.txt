4_PIC_operators_declaration -> Build the operators of the forward model and performs their autocalibration on the science data.

- 'executable.m' runs the operator declaration according to the parameters set in 'procedures/parameters.m'

- 'parameters.m' -> mainly set the parameters in the sections '%% Data' (to point to the spectral calibration) and '%% Shift determination' (to perform the autocalibration described in Appendix C via the variables 'flag_estimate_shift', 'exp_time', 'path_data' and 'path_calib_sensor', the other variables being generally untouched).
