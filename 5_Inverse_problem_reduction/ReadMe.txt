5_Inverse_problem_reduction -> Data reduction by inverse problem approach with the PIC forward model

- 'executable.m' -> run the reduction according to the parameters set in 'procedures/parameters.m'

- 'parameters.m' -> mainly set the parameters in the sections '%% Data', the rest being generally untouched (/!\ 'data.pos_center' is used to center the reduction on the star center given after the analysis of the reduction of the 'OBJECT_CENTER').

- 'analyze_centroid.m' -> determine the position of the star center in the reduction of the waffle mode acquisition (OBJECT_CENTER). Set the parameters in the section '%% Set parameters' (mainly the variables 'data.path', 'data.name' and 'list_lambda', the others being generally untouched. The display dynamics is set in the variable 'color_scale').