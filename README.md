# Eye_Tracking_Data_Processing_Pipeline

This selection of Matlab scripts reads asc files of SR research result files, detects fixation periods from markers in the file, summarizes samples to fixation locations, corrects data according to fixation events during the experiment, conducts ROIs analyses, calculates serveral eye tracking variables (e.g. number of fixations per trial, percent fixation time in ROI, average distance to ROI etc.) and visualizes data as heat maps.


ET_Data_Processing_Clear_Syllables.m
of SR research result files, detects fixation periods from markers in the file, summarizes samples to fixation locations

ET_Data_Processing_Clear_Syllables_Center_Fixes.m
If presentation script includes occasional fixation control events (e.g. in the center of the screen) this scropt reads out the fixations of the control fixations.

Correct_Data.m
Corrects experimental data according according to fixation events during the experiment

ROI_Processing_Clear_Syllables.m
Conducts ROIs analyses with previously specified ROIs, calculates serveral eye tracking variables (e.g.  percent fixation time in ROI, average distance to ROI etc.) and summarizes results per trial, run and subject

Heat_Map_tot_per_sub.m
Visualization of all trials (or separately over several experimental conditions) per subject as heat map

Heat_Map_tot_over_sub.m
Visualization of all trials (or separately over several experimental conditions) for the whole experiment as heat map
