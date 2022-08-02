#---#
clear all
close all

% addpath("../JCR/cruise_specific_functions")
addpath(strcat(pwd, "/cruise_specific_functions")) %FICE 2022: assumes code is run from ../HSAS_Processing 
addpath(strcat(pwd, "/rad_functions/"))
addpath(strcat(pwd, "/rad_functions/intwv"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa"))


input_parameters_hsas;
disp(DIN_HSAS)
