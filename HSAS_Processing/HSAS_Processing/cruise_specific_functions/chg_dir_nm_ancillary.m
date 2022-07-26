# change name of ancillary subdirectories	

clear all

addpath("~/Dropbox/Octave_functions/")
addpath("../hsas.source")
addpath("../hsas.source/cruise_specific_functions")
addpath("../hsas.source/rad_functions/")
addpath("../hsas.source/rad_functions/intwv")
addpath("../hsas.source/rad_functions/DISTRIB_fQ_with_Raman")
addpath("../hsas.source/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa")



input_parameters_hsas();

# rd original dir names
dirs = glob([DIR_GPS "2015*"])

for id = 1:1#2:length(dirs)

	# read original date
		o_dt = strsplit(dirs{id}, '/'){end};
		doy = str2num(o_dt(end-2:end));
		
		n_dt = datestr(datenum([str2num(o_dt(1:end-3)) 1 doy]), "yyyymmdd");

		new_dir = strrep(dirs{id}, o_dt, n_dt);
		
		system(["mv " dirs{id}, " ", new_dir]);
endfor


