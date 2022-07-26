# this is a script that contains all (I hope) the variables parameters for each cruise

#Turn off/on verbose reporting
	VBS = true;

# year of dataset
	DEF_YEAR = 2022;
	
# dates to be processed
	%DAY_START = "195"; %14 july %AMT proc uses jday
	%DAY_STOP = "202"; %21 july
	DAY_START = "20220713"; % 13 july % fice proc uses datestr
	DAY_STOP = "20220721";  % 21 july

# cruise name
	CRUISE = "FICE22";

# HSAS instruments serial numbers
	INSTRUMENT = "hsas";
	radiometers = {"ES", "LI", "LT"};#  (similar instrument must be listed one after the other)
	sn = {"2027", "2054", "464"};
	file_ext = {"H[ES][DE]", "*H[LS][DL]", "*H[LS][DL]", "SAT"}; # wildcards to read files for each instrument
	cal_files_pre = {"SAT2027.cal", "SAT2054.cal", "SAT464.cal"};
	cal_files_post = {"SAT2027.cal", "SAT2054.cal", "SAT464.cal"};
	
# Set wavelength range (note: this is not the raw wavelength - it is interpolated)
   wv = [350:2:860]';  #Consistent with JCR format

# main dir
MAIN_DIR = "/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/";   

% basename for calibrated files 
%FNBASE	= [num2str(DEF_YEAR) '-' doy] % filename base for AMT crusies: doy variable in proc1
FNBASE	= ['FRM4SOC2_FICE22_PML_']; %

### INSTRUMENT serial numbers for trios (at least) are hardcoded
### SINGLE HARDCODED PATH REMAINS TO THS DATA BELOW AND DEF_YEAR IS HARDCODED ###
OSX = 0;

# main directories
% DOUT_SUFFIX = "./Processed_final/"; # this is the supphix that is appended to the DATA_PATH dir to define the directory for the output processed files


DIR_GPS = [MAIN_DIR "xls/"];  # NOTE: each day of ancillary data must be separately stored in a directory with name yyyymmdd (e.g., 20150926)
GLOB_GPS = "/*fice22_metadata.xlsx";  % lat = 45.31425, lon = 12.508317

% case of AMT cruise
% DIR_ATT = [MAIN_DIR "Ship_uway/"]; # pitch and roll # NOTE: each day of ancillary data must be separately stored in a directory with name yyyymmdd (e.g., 20150926)
% GLOB_ATT = "/tsshrp.ACO"; 
% case of FICE22 - notheing needed;

DIR_WIND = [MAIN_DIR "xls/"];
GLOB_WIND = "/*meteo.xlsx";# glob pattern for wind data to concatenate after DATESTR

%DIR_SURF = [MAIN_DIR "xls/"]; % temp and humidity stored in wind struct for FICE
%GLOB_SURF = "/*meteo.xlsx";# glob pattern for surf data to concatenate after DATESTR
%DIR_SURF = [];


DIR_TEMPCORR = [MAIN_DIR "HyperSAS_config/Temperature/"];
FN_TEMPCORR_ES = "PML_8-010-20-thermal-0258.csv";# File name of temperature correction factors and their uncertainties
FN_TEMPCORR_LT = "PML_8-010-20-thermal-0223.csv";
FN_TEMPCORR_LI = "PML_8-010-20-thermal-0222.csv";

DIR_SLCORR = [MAIN_DIR "HyperSAS_config/Straylight/"];
FN_SLCORR_ES = "cp_straylight_SAT2027.txt";# File names of Straylight correction factors for FICE
FN_SLCORR_LT = "cp_straylight_SAT0464.txt";
FN_SLCORR_LI = "cp_straylight_SAT2054.txt";


DIN_HSAS = [MAIN_DIR "Raw/RawExtracted/"];
DOUT_HSAS_SUB = "Processed_TPcorrection/";
DOUT_HSAS = [MAIN_DIR DOUT_HSAS_SUB];

DIR_CAL = [MAIN_DIR "HyperSAS_config/"];
DIN_CALS_PRE = [DIR_CAL "Pre/"];
DIN_CALS_POST = [DIR_CAL "Post/"];

DIN_StrayLight = [MAIN_DIR "HyperSAS_config/Straylight/"];
DIN_Non_Linearity = [MAIN_DIR "HyperSAS_config/Non-linearity/"]; 

NL_files_pre = {"Pre/cp_radcal_SAT2027.txt", "Pre/cp_radcal_SAT2054.txt", "Pre/cp_radcal_SAT0464.txt"};
NL_files_post = {"Post/cp_radcal_SAT2027.txt", "Post/cp_radcal_SAT2054.txt", "Post/cp_radcal_SAT0464.txt"};

# -Flags to check if do non-linearity, temperature, Straylight corrections
FLAG_NON_LINEARITY = 0; % 0: do not apply correction, 1: apply correction
FLAG_TEMPERATURE = 0;
FLAG_STRAY_LIGHT = 0;

#-ACS data Path
FN_ACS = [MAIN_DIR, "AC-9/ACStoHSAS_sentine.txt"];

# Define names of functions needed to read GPS, HDG, PITCH, ROLL, TILT
FNC_RD_ATT = @rd_tsshrp_FICE;# function to read pitch and roll, heave
FNC_RD_GPS = @rd_gps_FICE;# function to read lat, lon, hdg, cog_deg, sog _m2s
FNC_RD_WIND = @rd_anemometer_FICE; # function to read wdir, wspd, airt, humid  % modify for xls
FNC_RD_SURF = @rd_oceanlogger_FICE; # function to read other met and surface data collecetd by the ship 


##########################################
## Parameters for L2 processing
# type of filtering applied to data
##########################################
%FILTERING  = 'continuous'; # 	L1_f = hsas_filter_sensors_using_Lt_data_v2(L1_f, L1_f, 'vaa_ths'); 
 FILTERING  = 'lowest'; # 		L1_f = hsas_filter_sensors_using_Lt_data(L1_f, L1_f, 25, 'vaa_ths')
 %FILTERING  = 'both'; # 		L1_f = hsas_filter_sensors_using_Lt_data_v2(L1_f, L1_f, 'vaa_ths');
		  			   # 		L1_f = hsas_filter_sensors_using_Lt_data(L1_f, L1_f, 25, 'vaa_ths');

#
# directory where the ship's underway data are stored
% DIR_SHIP = "/data/datasets/cruise_data/active/AMT24/Ship_data/Compress/Compress/days/";
#
# base dir fir L1 files
DIN_L1 = [MAIN_DIR DOUT_HSAS_SUB "L1/"];
#
# maximum tilt accepted 
MAX_TILT = 5; # [degrees]
#
# nominal viewing angle of Li sensor wrt zenith (depends on installation)
LI_VA = 40; # [degrees]






#----------------------------------
### Read THS data ###
%
% if OSX==1;
%    din_ths = ["/Volumes/Rivendell/AMT/raw_data/Hsas/SatCon_Extracted_oldcal/" num2str(doy) "/"];
% else
%    din_ths = [DATA_PATH "Hsas/Processed_final/Extracted/" num2str(doy) "/"];
% endif
%
% fnths   = glob( [din_ths "*SATTHS*dat"] );
% ths     = hsas_rd_ths(fnths);






% # above-water TRIOS instrument serial numbers
% 	sn = {"82C1", "ES"};
%     sn = {"8313", "LI"};
%     sn = {"8346", "LT"};
%
%
% # in-water TRIOS instrument serial numbers
%     sn = {"8508", "LU"};
%
% 	# Apply immersion factor for calibration in water vs air - currently uses factor calculated from default files provided by Marco
%     if_out = csvread('/data/datasets/cruise_data/active/AMT26/PML_optics/HEK_processing/code/inwater_cals/immersion_factors.csv');
	

# Method for interpolating to the same time step
TIME_INT_METHOD = "linear";

# Tilt filter
MAX_TILT_ACCEPTED_L1 = 5;

# SZA filter
MAX_SZA_L2 = 80; # degrees]
MIN_SZA_L2 = 10; # degrees]


# PHI filter
MAX_PHI_L2 = 170; # degrees]
MIN_PHI_L2 =  50; # degrees]
 



