# July 2016: Original code by gdal 
#
# 22nd March: 2018 Hayley Evers-King 
# Code to run L0 and L1 processing for above water hypersas and trios
#
# Spring 2019: Modified by gdal to process AMT24 data
#
# Oct 2019: New version to clean up and generalise by gdal
# Code set to process one file or one day at a time? ()


### Set up ###
clear all


warning('off');


#pkg load signal
pkg load netcdf

#struct_levels_to_print (0);
#debug_on_error(1);

warning off #Turn off warnings

addpath(strcat(pwd, "/cruise_specific_functions")) %assumes code is run from ../HSAS_Processing
addpath(strcat(pwd, "/rad_functions/"))
addpath(strcat(pwd, "/rad_functions/intwv"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa"))


# read input parameters for this cruise
input_parameters_hsas;


# Get arguments passed to function: INSTRUMENT switch is first argument xargs comes after
%fnin = argv; # tj - THIS NEEDS TO BE UNCOMMENTED 
fnin = {"20220714"}; % tj - TEMPORARY HARD CODING - 1st day of FICE2022
% fnin = {"20150916"};
% fnin = {"hsas", ...
%   		"20191017", ...
%   	  	"v1"			};#, ...
#  "/Volumes/Public/DY110_Public/Optics_group/Data/Underway_Rrs/20191025/20191025_Satcon_extracted_raw/2019-298-120943-HSL223I.dat"};

#fnin{1};
DATESTR = fnin{1};
VERSION   = "v1";#fnin{3};
 
doy = num2str(jday(datenum(DATESTR, "yyyymmdd")));
%doy = fnin{1}
 


% create strings with year, month and day
%dout   	= [DOUT_HSAS DATESTR "/" ];
dout_L0   	= [DOUT_HSAS "L0/" DATESTR "/" ];
dout_L1   	= [DOUT_HSAS "L1/" DATESTR "/" ];
% dbase  	= strsplit(fnin{4}, {"/"}){end};
% doy  	= strsplit(dbase, {"-"}){end-2};
%Y0     	= datenum(DEF_YEAR, 1, 1);
%[Y, M, D] = datevec(datenum(DATESTR, "yyyymmdd"));

fnbase 	= [num2str(DEF_YEAR) '-' doy];
YYYY = sprintf("%4u", DEF_YEAR);
[M,D] = jday2mmdd(DEF_YEAR,str2num(doy)); % CODE currently breaks here
MM = sprintf("%02u", M);
DD = sprintf("%02u", D);
DATESTR2 = [num2str(DEF_YEAR) doy];



 
# Set up directories and file names
if strcmp(INSTRUMENT,'hsas')
   din    = [DOUT_HSAS "Calibrated/" DATESTR "/" ];      # <<<<<<<<<<< CHECK THIS
   disp('Processing hypersas....');
   fflush(stdout);
   
   in_tag 	= 'Hsas';	
   outbase = [YYYY '_' MM '_' DD '_' in_tag];
 
%
% elseif strcmp(INSTRUMENT,'triosAW')
%    inst_type = "TRIOS";
%    din    = fnin{4}(1:strfind(DATA_PATH,'/')(end));
%    disp('Processing trios....')
%    dout = [DATA_PATH DOUT_SUFFIX];
%    dbase  = strsplit(fnin{4}, {"/"}){end};# <<<<<<<<<<< CHECK THIS
%    dbase  = strsplit(dbase, {"_"}){end-2};# <<<<<<<<<<< CHECK THIS
%    [Y,M,D]= datevec(dbase, "yyyy-mm-dd");
%    Y0     = datenum(Y,1,1);
%    doy    = datenum(Y,M,D) - Y0 + 1;
%    fnbase = strsplit(fnin{4}, {"/","."}){end-1};# <<<<<<<<<<< CHECK THIS
%    in_tag = 'Trios';
%    inst_name = strsplit(fnbase,"_"){2};# <<<<<<<<<<< CHECK THIS
%
% elseif strcmp(INSTRUMENT,'triosIW')
%    inst_type = "TRIOS";
%    din    = fnin{4}(1:strfind(DATA_PATH,'/')(end)); # <<<<<<<<<<< CHECK THIS
%    disp('Processing trios....')
%    dout = [DATA_PATH DOUT_SUFFIX];
%    dbase  = strsplit(fnin{4}, {"/"}){end};# <<<<<<<<<<< CHECK THIS
%    dbase  = strsplit(dbase, {"_"}){end-2};# <<<<<<<<<<< CHECK THIS
%    [Y,M,D]= datevec(dbase, "yyyy-mm-dd");
%    Y0     = datenum(Y,1,1);
%    doy    = datenum(Y,M,D) - Y0 + 1;
%    fnbase = strsplit(fnin{4}, {"/","."}){end-1};# <<<<<<<<<<< CHECK THIS
%    in_tag = 'Trios';
%    inst_name = strsplit(fnbase,"_"){2};# <<<<<<<<<<< CHECK THIS

else 
	keyboard

endif



#----------------------------------
### Read GPS data ###
#
 %here we should use functions name in input parameters to read gps data 
	fn_gps = glob([DIR_GPS, DATESTR, GLOB_GPS]){1}; % fill with lat lon from platform
	gps = FNC_RD_GPS(  fn_gps   );
	
	
#----------------------------------
### Read TILT AND ROLL data ###
	% fn_ths = glob([DIR_ATT, DATESTR2, GLOB_ATT]){1};
	% ths = FNC_RD_ATT(  fn_ths   );
	
	% ths.tilt = cmp_tilt(ths.roll, ths.pitch); # fill with tilt =0
	
	
#----------------------------------
### Read WIND data ###
# here we should use functions name in input parameters to read gps data 
	fn_wind = glob([DIR_WIND, DATESTR, GLOB_WIND]){1}; % file name of wind xlsx 
	wind = FNC_RD_WIND(  fn_wind   ); % wind data are not in utc - check bias correction in xls file
	keyboard
#----------------------------------
### Read other met and surface data ###
# here we should use functions name in input parameters to read gps data 
if exist(DIR_SURF)
	fn_surf = glob([DIR_SURF, DATESTR2, GLOB_SURF]){1};
	surf = FNC_RD_SURF(  fn_surf   ); % similar to wind
endif


#---------------------------------
####### Read non-linearity correction coefficients ##########
#---Moved to Calibration step
pkg load io
#---radiometer related to sn
#rad_sn = cell2struct(sn,radiometers,2);
#coeff_LI = xlsread(DIN_Non_Linearity,rad_sn.LI);
#coeff_LT = xlsread(DIN_Non_Linearity,rad_sn.LT);
#coeff_ES = xlsread(DIN_Non_Linearity,rad_sn.ES);
#non_linearity_coeff = struct('coeff_LI',coeff_LI(:,1:2),'coeff_LT',coeff_LT(:,1:2),'coeff_ES',coeff_ES(:,1:2));



############################################################
###### HSAS processing #####################################

if strcmp(INSTRUMENT, 'hsas')
	
	
	iinstr = 1;	# this is a counter for instrument with similar names (e.g., two LT radiometers)
	 
  
	for irad = 1:length(radiometers)
			
       # intialize structure		
		clear rad	
			
	   # create structure with SN and radiometer type (ES/LI/LT)
		sn_rad = { sn{irad}, radiometers{irad} };
		
 	   # read raw data
		% fn = glob( [din fnbase "*" file_ext{irad} sn{irad} "*dat"] );
		fn = glob( [din fnbase "*-[0-2]*" file_ext{irad} sn{irad} "*dat"] );  ####################################################
		rad.raw = hsas_rd_many(sn_rad, fn, VBS); 

    
 	   # correct for dark counts
 	   	rad.raw_nodk = correct_dk(sn_rad, rad.raw, VBS); 
                   
     # non-linearity correction
      #rad = correct_non_linearity(rad_sn, rad.raw.sn,non_linearity_coeff,rad);

 	   # interpolate using wv as a new wavelength set
 	   	rad.raw_intwv = intwv(sn_rad, rad.raw_nodk, wv, VBS);
		

 	   # check if the previous type of instrument has already been processed (requires that similar instrument are listed one after the other in input_parameters.m)
 	   if irad>1 & strcmp(radiometers{irad}, radiometers{irad-1})
 		   iinstr = iinstr + 1; # and if it has, increments the counter
 	   else
 		   iinstr = 1;
 	   endif

	   # create structure for this instrument
	   	if iinstr == 1
	   		eval([radiometers{irad} " = rad;"]);
		else
			eval([radiometers{irad} num2str(iinstr) " = rad;"]);
		endif
    
	endfor
	

   
endif



############################################################
###### Above-water TRIOS processing ########################

	if strcmp(INSTRUMENT, 'triosAW')
	   % sn = {"82C1", "ES"};

	   # read raw ES data
	   fn = glob( [din fnbase "*H[ES][DE]" sn{1} "*dat"] );# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	   ES.raw = rd_trios(fnin{4});# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	   ES.raw.sn = sn{1};# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS

	   # correct for dark counts
	   #ES.raw_nodk = correct_dk(sn, ES.raw, VBS);

	   # interpolate using wv as a new wavelength set
	   ES.raw_intwv = intwv(sn, ES.raw, wv, VBS);
	endif

	if strcmp(INSTRUMENT, 'triosAW')
	    % sn = {"8313", "LI"};

	    # read raw LI data
	    if VBS, disp('reading LI file'); fflush(stdout); endif
	    LI.raw = rd_trios(fnin{5});# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	    LI.raw.sn = sn{1};# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS

	    # correct for dark counts
	    #LI.raw_nodk = correct_dk(sn, LI.raw, VBS); # previously commented - why? hsas only

	    # interpolate using wv as a new wavelength set
	    LI.raw_intwv = intwv(sn, LI.raw, wv, VBS);
	endif

	if strcmp(INSTRUMENT,'triosAW')
	    % sn = {"8346", "LT"};

	    # read raw Lt data
	    LT.raw = rd_trios(fnin{6}); # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	    LT.raw.sn = sn{1};

	    # correct for dark counts
	    #LT.raw_nodk = correct_dk(sn, LT.raw, VBS); # previously commented - why? hsas only

	    # interpolate using wv as a new wavelength set  
	    LT.raw_intwv = intwv(sn, LT.raw, wv, VBS);
	endif

############################################################
###### Below-water TRIOS processing ########################
	#----------------------------------
	### Read Lu ###
	if strcmp(INSTRUMENT,'triosIW')
	    % sn = {"8508", "LU"};

	    # read raw Lu data
	    LU.raw = rd_trios(fnin{4});   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	    LU.raw.sn = sn{1};

	    # correct for dark counts
	    #LU.raw_nodk = correct_dk(sn, LU.raw, VBS); # previously commented - why? hsas only

	    # interpolate using wv as a new wavelength set  
	    LU.raw_intwv = intwv(sn, LU.raw, wv, VBS);
    
	    % # Apply immersion factor for calibration in water vs air - currently uses factor calculated from default files provided by Marco
	    % if_out=csvread('/data/datasets/cruise_data/active/AMT26/PML_optics/HEK_processing/code/inwater_cals/immersion_factors.csv');

	    # interpolate using wv as a new wavelength set  
	    if_intwv = interp1(LU.raw.wv, if_out(13:162,3)', wv);# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS
	    LU.raw_intwv.data = LU.raw_intwv.data.*if_intwv';
	endif










#----------------------------------
### Create L0 structure ###

L0.wv = wv;
if strcmp(INSTRUMENT,'triosIW')
   L0.time = LU.raw_intwv.time;
else
   L0.time = LT.raw_intwv.time;
endif

#----------------------------------
# interpolate all data to LT time step

if strcmp(INSTRUMENT,'triosIW')
   # IN WATER RADIOMETER ONLY MEASURES LU
   L0.instr.Lu = hsas_int2Lt(LU.raw_intwv, LU.raw_intwv.time);
   L0.instr.Lu.roll = interp1(ths.time, ths.roll, L0.time, TIME_INT_METHOD);
   L0.instr.Lu.pitch = interp1(ths.time, ths.pitch, L0.time, TIME_INT_METHOD);
   L0.instr.Lu.tilt = interp1(ths.time, ths.tilt, L0.time, TIME_INT_METHOD);
   #L0.instr.Lu.vaa_ths = interp1(ths.time, ths.compass, L0.time, TIME_INT_METHOD);
   try 
      L0.instr.Lu.vaa_ths = interp1(ths.time, ths.compass, L0.time, TIME_INT_METHOD);
   catch
      L0.instr.Lu.vaa_ths = interp1(gps.hdg_time, gps.hdg, L0.time, TIME_INT_METHOD);
   end_try_catch

   L0.instr.Lu.sn = LU.raw.sn;
   L0.tilt = abs([L0.instr.Lu.tilt]);
   L0.roll = abs([L0.instr.Lu.roll]);
   L0.pitch = abs([L0.instr.Lu.pitch]);
   
endif

% if ~strcmp(INSTRUMENT,'triosIW')
if strcmp(INSTRUMENT,'hsas')
	
   # ABOVE WATER RADIOMETER MEASURES LT,LI,ES
   # Es
   L0.instr.Es = hsas_int2Lt(ES.raw_intwv, LT.raw_intwv.time);
   L0.instr.Es.roll = interp1(ths.time, ths.roll, L0.time, TIME_INT_METHOD);
   L0.instr.Es.pitch = interp1(ths.time, ths.pitch, L0.time, TIME_INT_METHOD);
   L0.instr.Es.tilt = interp1(ths.time, ths.tilt, L0.time, TIME_INT_METHOD);
   L0.instr.Es.temp_pcb = interp1(ES.raw.time, ES.raw.TEMP_degC, LT.raw_intwv.time);
   % try
   %    L0.instr.Es.vaa_ths = interp1(ths.time, ths.compass, L0.time, TIME_INT_METHOD);
   % catch
      L0.instr.Es.vaa_ths = interp1(gps.time, gps.hdg, L0.time, TIME_INT_METHOD);
   % end_try_catch

   L0.instr.Es.sn = ES.raw.sn;


   # Li
   L0.instr.Li = hsas_int2Lt(LI.raw_intwv, LT.raw_intwv.time);
   L0.instr.Li.roll = interp1(ths.time, ths.roll, L0.time, TIME_INT_METHOD);
   L0.instr.Li.pitch = interp1(ths.time, ths.pitch, L0.time, TIME_INT_METHOD);
   L0.instr.Li.tilt = interp1(ths.time, ths.tilt, L0.time, TIME_INT_METHOD);
   L0.instr.Li.temp_pcb = interp1(LI.raw.time, LI.raw.TEMP_degC, LT.raw_intwv.time);
   % try
 %      L0.instr.Li.vaa_ths = interp1(ths.time, ths.compass, L0.time, TIME_INT_METHOD);
 %   catch
      L0.instr.Li.vaa_ths = interp1(gps.time, gps.hdg, L0.time, TIME_INT_METHOD);
   % end_try_catch

   L0.instr.Li.sn = LI.raw.sn;


   # Lt
   L0.instr.Lt = hsas_int2Lt(LT.raw_intwv, LT.raw_intwv.time);
   L0.instr.Lt.roll = interp1(ths.time, ths.roll, L0.time, TIME_INT_METHOD);
   L0.instr.Lt.pitch = interp1(ths.time, ths.pitch, L0.time, TIME_INT_METHOD);
   L0.instr.Lt.tilt = interp1(ths.time, ths.tilt, L0.time, TIME_INT_METHOD);
   L0.instr.Lt.temp_pcb = interp1(LT.raw.time, LT.raw.TEMP_degC, LT.raw_intwv.time);
   % try
   %    L0.instr.Lt.vaa_ths = interp1(ths.time, ths.compass, L0.time, TIME_INT_METHOD);
   % catch
      L0.instr.Lt.vaa_ths = interp1(gps.time, gps.hdg, L0.time, TIME_INT_METHOD);
   % end_try_catch

   L0.instr.Lt.sn = LT.raw.sn;


   # Calculate tilt/pitch/roll filter
   L0.tilt = max( abs( [L0.instr.Es.tilt, L0.instr.Li.tilt, L0.instr.Lt.tilt,  L0.instr.Es.tilt - L0.instr.Li.tilt, L0.instr.Es.tilt - L0.instr.Lt.tilt, L0.instr.Li.tilt - L0.instr.Lt.tilt]' )   )';
   L0.roll = max( abs( [L0.instr.Es.roll, L0.instr.Li.roll, L0.instr.Lt.roll,  L0.instr.Es.roll - L0.instr.Li.roll, L0.instr.Es.roll - L0.instr.Lt.roll, L0.instr.Li.roll - L0.instr.Lt.roll]' )   )';
   L0.pitch = max( abs( [L0.instr.Es.pitch, L0.instr.Li.pitch, L0.instr.Lt.pitch,  L0.instr.Es.pitch - L0.instr.Li.pitch, L0.instr.Es.pitch - L0.instr.Lt.pitch, L0.instr.Li.pitch - L0.instr.Lt.pitch]' )   )';



   # add the apparent wind data

	L0.appwind_spd = interp1(wind.time, wind.ws, L0.time, TIME_INT_METHOD); # [m/s]
	L0.appwind_dir = interp1(wind.time, wind.wdir, L0.time, TIME_INT_METHOD); # [degrees]
	
	if isvarname('surf')
		fld = fieldnames(surf);
		for ifld = 1:length(fld)
			L0.surf_met.(fld{ifld}) = interp1(surf.time, surf.(fld{ifld}), L0.time, TIME_INT_METHOD); # [degC]
		endfor
	else
		if isfield(wind, 'airt')
			L0.surf_met.airtemp = interp1(wind.time, wind.airt, L0.time, TIME_INT_METHOD); # [degC]
			L0.surf_met.humid = interp1(wind.time, wind.humid, L0.time, TIME_INT_METHOD); # [%] relative humidity
		endif
	endif
	
	
	
	
endif




#----------------------------------
# Header

L0.hdr.site_code = "A";

for i_instr = 1:length(fieldnames(L0.instr))
   L0.hdr.method_code(i_instr,:) = "04A020";
endfor

if strcmp(INSTRUMENT,'hsas')
   L0.hdr.INSTRUMENT_id = "HyperOCR";
else
   L0.hdr.INSTRUMENT_id = "Trios";
endif

L0.hdr.station_ref = CRUISE;

#### Straylight and temperature correction ####
#---Stray Light Correction (Moved to Calibration step)
# L0 = hsas_straylight_correct(L0);

#---Temperature Correction
if FLAG_TEMPERATURE == 1
L0 = hsas_temperature_correct(L0);
endif

### Create L1 structure ###

# Tilt filter
%
% max_tilt_accepted = 50;

L1.hdr.max_tilt_accepted = MAX_TILT_ACCEPTED_L1;

if strcmp(INSTRUMENT,'triosIW')
   L1 = hsas_mk_L1_filtered(L0, MAX_TILT_ACCEPTED_L1,"IW");
else
   L1 = hsas_mk_L1_filtered(L0, MAX_TILT_ACCEPTED_L1,"AW"); 
endif


##### remove old lat lon fields
#    L1 = rmfield(L1, {"lat","lon"});


##### add gps data

L1.gps.lat = nan(size(L1.time));
L1.gps.lon = nan(size(L1.time));
L1.gps.time = nan(size(L1.time));
L1.gps.hdg = nan(size(L1.time));
L1.gps.sog_m2s = nan(size(L1.time)); # speed over ground [m/s]
L1.gps.cog_deg = nan(size(L1.time)); # course over ground [degrees]

L1.gps.lat = interp1(gps.time, gps.lat, L1.time);
L1.gps.lon = interp1(gps.time, gps.lon, L1.time);
L1.gps.hdg = interp1(gps.time, gps.hdg, L1.time);
L1.gps.sog_m2s = interp1(gps.time, gps.sog_m2s, L1.time);
L1.gps.cog_deg = interp1(gps.time, gps.cog_deg, L1.time);

#----computer rho
L1 = hsas_cmp_rho_2pars(L1)

### Save files ###

# Binary format
%dout_L0 = [dout "L0/" ];
if ~exist(dout_L0)
    mkdir(dout_L0);
endif      
fnout_L0 = [dout_L0 outbase "_L0_" VERSION];
save("-v7", [fnout_L0 ".mat"],  "L0", "doy");



%dout_L1 = [dout "L1/" ];
if ~exist(dout_L1)
    mkdir(dout_L1);
endif
fnout_L1 = [dout_L1 outbase "_L1_" VERSION];
save("-v7", [fnout_L1 ".mat"],  "L1", "doy");



