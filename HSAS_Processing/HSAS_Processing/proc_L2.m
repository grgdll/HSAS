# Hayley Evers-King 8th Feb 2017
# Code to run L2 processing for above water hypersas and trios
# modified by gdal to process AMT24 data (22 Nov 2018)

### Set up ###
clear all

% global main_dir DIR_SHIP

struct_levels_to_print (0);

warning off #Turn off warnings

% addpath("~/Dropbox/Octave_functions/")
addpath(strcat(pwd, "/cruise_specific_functions")); % assumes code is run from ../HSAS_Processing
addpath(strcat(pwd, "/rad_functions/"));
addpath(strcat(pwd, "/rad_functions/intwv"));
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman"));
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa"));


# read input parameters for this cruise
input_parameters_hsas;

%fnin = argv; # tj - THIS NEEDS TO BE UNCOMMENTED 
fnin = {"20220714"}; % tj - TEMPORARY HARD CODING - 1st day of FICE2022
% fnin = {"20191019"};
# fnin ={"hsas",...
#   "/data/lazarev1/backup/cruise_data/AMT24/DallOlmo/HSAS/Satcon_extracted/Physical_units/",...
#   "v1",...
#   "continuous", ...
# #   "/data/lazarev1/backup/cruise_data/AMT24/DallOlmo/HSAS/Satcon_extracted/Physical_units/../../Processed_final/L1/275/2014_10_02_Hsas_L1_v1.mat"}   
#    "/data/lazarev1/backup/cruise_data/AMT24/DallOlmo/HSAS/Satcon_extracted/Physical_units/../../Processed_final/L1/272/2014_09_29_Hsas_L1_v1.mat"}

% instrument = fnin{1};
% data_path  = fnin{2}


disp('----------------------------------------------------------------------------------------------------------');
disp(fnin{1});
disp('----------------------------------------------------------------------------------------------------------');
fflush(stdout);


DATESTR = fnin{1};
VERSION   = "v1";
 
%doy = num2str(jday(datenum(DATESTR, "yyyymmdd")));
doy = fnin{1}

[DIN_L1 DATESTR "/*mat" ]
filename = glob([DIN_L1 DATESTR "/*mat" ]){1};
if isempty(filename)
    disp('WARNING: No input file exiting.')
    exit
endif



# for use in sub-functions
main_dir  = [];

### Set up directories and filenames ###

din = DIN_L1; #filename(1:strfind(fnin{2}, '/')(end));
fnbase = [DIN_L1 DATESTR "/"];

if strcmp(INSTRUMENT,'hsas');
   dout = strrep(fnbase, 'L1','L2');
   dout_nc = [dout 'netcdf/'];
endif


% if strcmp(fnin{1},'triosAW');
%    dout = [data_path '/TriosAW/Processed_final/L2/'];
% endif
% if strcmp(fnin{1},'triosIW');
%    dout = [data_path '/TriosIW/Processed_final/L2/'];
% endif

# create output dir if it does not exist
if ~exist(dout)
   system(["mkdir -p " dout]);
endif

if ~exist(dout_nc)
   system(["mkdir -p " dout_nc]);
endif



# prepare file names for output
fnout = strrep(filename, "L1", "L2");

# Check if already processed
if exist(fnout)
     disp('file already processed: exiting.')
     #exit
endif



### Load file ###
load(filename);

# FOR IW DATA, I NEED TO INTERP ES FROM AW DATA, DO THAT HERE
% if strcmp(INSTRUMENT,'triosIW');
%    fnin2=strrep(fnin{5},"TriosIW","Hsas");
%    fnin2=strrep(fnin2,"Trios","Hsas");
%    disp(fnin2)
%    IW_TIME   = L1.instr.Lu.time;
%    IW_WAVE   = L1.wv;
%    AW_STRUCT = load(fnin2);
%    AW_TIME   = AW_STRUCT.L1.instr.Es.time;
%    AW_DATA   = AW_STRUCT.L1.instr.Es.data;
%    AW_WAVE   = AW_STRUCT.L1.wv;
%    [AW_WAV_2D,AW_TIME_2D] = meshgrid(AW_WAVE,AW_TIME);
%    [IW_WAV_2D,IW_TIME_2D] = meshgrid(IW_WAVE,IW_TIME);
%    IW_DATA   = interp2(AW_WAV_2D,AW_TIME_2D,AW_DATA,IW_WAV_2D,IW_TIME_2D);
%
%
%    # make copy of Lu to hold Es data
%    L1.instr.Es = L1.instr.Lu;
%    L1.instr.Es.data = IW_DATA;
% endif









### compute average Li, Lt, Es ###

flds = fieldnames(L1.instr);

### apply filter for tilt ###

max_tilt_accepted = MAX_TILT; 
L1_f.hdr.max_tilt_accepted  = max_tilt_accepted;
L1_f.conf.max_tilt_accepted = max_tilt_accepted;  # yes this is redundant, but I did it to make sure all ancillary info is stored in L2.conf

if strcmp(INSTRUMENT,'triosIW')
   L1_f = hsas_mk_L1_filtered(L1, max_tilt_accepted, "IW");
else
   L1_f = hsas_mk_L1_filtered(L1, max_tilt_accepted, "AW");
endif

### Add fields to L1_f structure ###
L1_f.doy              = num2str(doy);
L1_f.files.L1         = filename;
L1_f.files.input_argv = {INSTRUMENT, DIN_L1, FILTERING};



####







keyboard 
### Create L2 structure ###
L2 = L1_f;

if strcmp(FILTERING,'lowest');
   # filters to leave only lowest n%
   if strcmp(INSTRUMENT,'triosIW')
      L1_f     = hsas_filter_sensors_using_Lu_data(L1_f, L1_f, 25, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lu_data(L2, L2, 25, 'vaa_ths');
   else
      L1_f     = hsas_filter_sensors_using_Lt_data(L1_f, L1_f, 25, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lt_data(L2, L2, 25, 'vaa_ths');
   endif
   FILT_TAG = '_Low';
   
elseif strcmp(FILTERING,'continuous');
   # this is the version for continuous data; temporal filtering
   if strcmp(INSTRUMENT,'triosIW')
      L1_f     = hsas_filter_sensors_using_Lu_data_v2(L1_f, L1_f, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lu_data_v2(L2, L2, 'vaa_ths');
   else
      L1_f     = hsas_filter_sensors_using_Lt_data_v2(L1_f, L1_f, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lt_data_v2(L2, L2, 'vaa_ths');
   endif
   FILT_TAG = '_Cont';
   
elseif strcmp(FILTERING,'both');
   if strcmp(INSTRUMENT,'triosIW')
      L1_f     = hsas_filter_sensors_using_Lu_data_v2(L1_f, L1_f, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lu_data_v2(L2, L2, 'vaa_ths');
      L1_f     = hsas_filter_sensors_using_Lu_data(L1_f, L1_f, 25, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lu_data(L2, L2, 25, 'vaa_ths');
   else
      L1_f     = hsas_filter_sensors_using_Lt_data_v2(L1_f, L1_f, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lt_data_v2(L2, L2, 'vaa_ths');
      L1_f     = hsas_filter_sensors_using_Lt_data(L1_f, L1_f, 25, 'vaa_ths');
      L2       = hsas_filter_sensors_using_Lt_data(L2, L2, 25, 'vaa_ths');
   endif
   
   FILT_TAG = '_both';
   
else
   FILT_TAG = '';
   
endif

keyboard 
% notes for begining work on 02/08
% (i) proc L2 currently runs as far as final nc write function (very nearly the end)
% (ii) hsas_filter_sensors_using_Lt_data result in a very large reduction in the size of L1_f and L2 (~1900 to 38 elements) We need to check why this is.
% (iii) There is underlying issue that time fields (L0, L1, L2) contain non-unqiue values: can be traced back to raw data
%

#### Compute geometries

# Compute sun position 
correction4PCtime = 0;

[solar_azimuth solar_elevation] = SolarAzEl(L1_f.instr.Es.time+correction4PCtime, L1_f.gps.lat, L1_f.gps.lon, Altitude=0);

L2.saa = L1_f.saa = solar_azimuth;
L2.sza = L1_f.sza = 90 - solar_elevation;

# Compute sensor position 
if strcmp(INSTRUMENT,'hsas')
   Li_angle = LI_VA;      
   L2.vza = L1_f.vza = Li_angle*ones(size(L1_f.gps.lat)) + L1_f.tilt;  # DO I NEED TO [ADD or SUBTRACT?] [tilt or pitch]?
   
   if isfield(L1_f, "gps")
      L2.vaa = L1_f.vaa = L1_f.gps.hdg;
   else 
      L2.vaa = L1_f.vaa = nan(size(L1_f.vza));
   endif
   
endif

% if strcmp(INSTRUMENT,'triosAW');
%     L2.vza = L1_f.vza = 40*ones(size(L2.gps.lat)) + L2.instr.Lt.tilt;  # DO I NEED TO [ADD or SUBTRACT?] [tilt or pitch]?
%     if isfield(L1_f, "gps")
%         L2.vaa = L1_f.vaa = L2.gps.hdg;
%     else
%         L2.vaa = L1_f.vaa = nan(size(L2.vza));
%     endif
% endif
%
% if strcmp(INSTRUMENT,'triosIW');
%     L2.vza = L1_f.vza = 0*ones(size(L2.gps.lat));
%     if isfield(L1_f, "gps")
%         L2.vaa = L1_f.vaa = L2.gps.hdg;
%     else
%         L2.vaa = L1_f.vaa = nan(size(L2.vza));
%     endif
% endif

# Compute relative azimuth between sensor and sun (phi)
if ~isfield(L1_f,'phi') % case of amt cruise
	L2.phi = L1_f.phi = L2.vaa - L2.saa;
	L2.phi(L2.phi>180) = L1_f.phi(L1_f.phi>180) = abs(360-abs(L2.phi(L2.phi>180))); # this is because the phi-axis of the fQ table ranges from 0 to 180.
	L2.phi(L2.phi<  0) = L1_f.phi(L1_f.phi<  0) = abs(L2.phi(L2.phi<0));
else % case of fice 2022
	L2.phi = L1_f.phi;
endif		        

### Read APPARENT????? wind speed and direction (this step is ship- and cruise-specific) 

   % L2 = rd_true_wind(L2);
   [~, L2.true_wind_dir, L2.true_wind_spd, nw, nwam, nwpm, nwf] = cmp_truew(sel=3, L2.gps.cog_deg, L2.gps.sog_m2s, L2.appwind_dir, zlr=0, L2.gps.hdg, L2.appwind_spd, wmis=nan(length(L2.gps.hdg),5) );

   L1_f.appwind_spd = L2.appwind_spd; # consistency
   L1_f.gps.sog_m2s = L2.gps.sog_m2s; # consistency


% ### Read atmospheric pressure and relative humidity    <<<<<<<<<<<<<<<<<   <<<<<<<   <<<<<<   NOT SURE WHY WE NEED THIS
%    fnocean = [DIR_SHIP 'oceanlogger/oceanlogger.' L2.doy];
%    tmpocean = rd_oceanloggerJCR(fnocean);
%
%    #check if there are significant amount of missing data
%       max_delta_seconds_ocean = max(diff(tmpocean.time)*1440*60);
%       if max_delta_seconds_ocean>60
%          disp(["WARNING: there is a LARGE GAP in the data from the oceanlogger." L2.doy " file: replace interpolation with intersection (or something else)"] );
%          fflush(stdout);
%          L2.WARNINGS.(['oceanlogger_' L2.doy]) = ["WARNING: there is a LARGE GAP in the data from the oceanlogger." L2.doy " file: replace interpolation with intersection (or something else)"];
%
%       endif
%
%     # interpolate ocean data to radiometric data
%         flds = fieldnames(tmpocean);
%         for ifld = 1:length(flds)
%             if strcmp(flds{ifld}, "time")
%                 continue
%             endif
%             L2.oceanlogger.(flds{ifld}) = interp1(tmpocean.time, tmpocean.(flds{ifld}), L2.time);
%         endfor
%


################################################################################
### Quality control ### NO MORE L1_F!

# SZA<=80
ikeep = find(L1_f.sza<=MAX_SZA_L2 & L1_f.sza>=MIN_SZA_L2);

if ikeep != 0.0
   L2 = hsas_filter(L2, ikeep);
   
else
   # Save file
   disp('-----------------------------No valid SZA, saving.')
   fnout = strrep(fnout, ".mat", [FILT_TAG ".mat"]);
   save("-v7", fnout, "L2", "L1_f");
   disp(["Written file " fnout FILT_TAG  ".mat"]);
   fflush(stdout);
  % exit
   
endif


#80<PHI>170
ikeep = find(L2.phi<MAX_PHI_L2 & L2.phi>MIN_PHI_L2);

if ikeep != 0.0
   L2 = hsas_filter(L2, ikeep); 
   
else
   # Save file
   disp('No valid PHI, saving.')
   save("-v7", fnout, "L2", "L1_f");
   disp(["Written file " fnout FILT_TAG  ".mat"]);
   fflush(stdout);
   % exit
   
endif



## Lt(750)>4 mW/(m2 nm sr)??? # commented in original code
# ikeep = find(L2.instr.Lt.data(:,find(L2.wv>=750,1))<1.5);
# L2 = hsas_filter(L2, ikeep);

### Calculate rho so that Lw in the NIR is zero and flat ###

if strcmp(INSTRUMENT,'triosIW');
   disp("No rho correction");
else
 %  L2 = hsas_cmp_rho_2pars(L2); % open ocean rho method
    L2 = hsas_extract_rho_AeronetOC(L2, 'appwind_spd'); % AAOT (tower) - apparent windspeed can be used
endif
	

# ### Calculate Lw = Lt - rho*Li
# if strcmp(fnin{1},'triosIW');
#    disp("Calculating Lw")
# 
#    ## Perform self shading correction as per Mueller et al., 2003 - ocean optics protocol, with some assumptions, including b<<a ##
# 
#    # Read in, sum  (ap plus mean acdom), and interpolate/extract absorption data from underway system
#    abs_file = '/data/datasets/cruise_data/active/AMT26/PML_optics/HEK_processing/data/TriosIW/ss_corr/abs_uw.oct';
#    load(abs_file);
#    
#    ap = cdom.hr.ap;
#    ac = (cdom.hr.ay+cdom.hr.cy)/2;
#    z1=zeros(size(ap,1),25);
#    z2=zeros(size(ap,1),55); 
#    ap = [ z1 ap z2];
#    ac = [ z1 ac z2];
#    aw_out=csvread('water_abs.csv');
#    aw_wv=[350:5:860]';
#    aw=interp1(aw_wv,aw_out(:,2),L2.wv(1,:));
# 
#    tot_a = ap + ac + aw;
#    tot_a(tot_a==0)=NaN;
# 
#    uw_wv = [350:2:860];
#    iw_wv = L2.wv; 
#    uw_time = cdom.hr.time;
#    iw_time = L2.time;
#    [iw_time_2D,iw_wv_2D] = meshgrid(iw_time,iw_wv);
#    [uw_time_2D,uw_wv_2D] = meshgrid(uw_time,uw_wv);
#    int_a = interp2(uw_time_2D, uw_wv_2D, tot_a', iw_time_2D, iw_wv_2D)';
# 
#    # Read in underway ship measurements of relative humidity and air pressure needed as input to the Gregg and Carder model used below (wind speed comes from processed file already)
#    uway_ship_file = '/users/modellers/gdal/cruise_data/active/AMT26/PML_optics/source/source/uway/output/Final/amt26_underway_optics_v7.mat';
#    load(uway_ship_file);
#    uws_time = amt26.uway.time;
#    uws_pres = (amt26.uway.baro1+amt26.uway.baro2)/2;
#    uws_relh = (amt26.uway.humidity1+amt26.uway.humidity2)/2;
#    uws_ws   = L2.ws;
#    int_pres = interp1(uws_time(isfinite(uws_pres)), uws_pres(isfinite(uws_pres)), iw_time);
#    int_relh = interp1(uws_time(isfinite(uws_relh)), uws_relh(isfinite(uws_relh)), iw_time);
# 
#    for i_Lu = 1:size(L2.instr.Lu.data,1)
# 
#    # Calculate likely clear sky Es fromDATESTR Gregg and Carder, plus ratio of direct:diffuse
#       [totI(i_Lu,:), dirI(i_Lu,:), difI(i_Lu,:)] = radmod_adap_hek(L2.time(i_Lu), L2.gps.lat(i_Lu), L2.gps.lon(i_Lu), L2.ws(i_Lu), int_pres(i_Lu), int_relh(i_Lu));
# 
#    # Compare clear sky Es with measurement Es and get scale factor
#       diff_tot(i_Lu,:)=L2.instr.Es.data(i_Lu,:)./totI(i_Lu,:);
#   
#    # Calculate h (scaled ratio of Esky (diffuse)/Esun (direct))
#       sf_h(i_Lu,:)=(difI(i_Lu,:)./dirI(i_Lu,:)).*diff_tot(i_Lu,:);
# 
#    # Compute g - ratio of sensor: instrument diameter - provided by Marco T 
#       g = 0.05/0.25;
# 
#    # Calculate K'sun
#       if (L1_f.sza(i_Lu) >=10) && (L1_f.sza(i_Lu) <20);
#          ksun0 = 3.14; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.56; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu)>=20) && (L1_f.sza(i_Lu) <30);
#          ksun0 = 3.05; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.49; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu) >=30) && (L1_f.sza(i_Lu) <40);
#          ksun0 = 2.94; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.39; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu) >=40) && (L1_f.sza(i_Lu) <50);
#          ksun0 = 2.80; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.28; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu) >=50) && (L1_f.sza(i_Lu) <60);
#          ksun0 = 2.64; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.15; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu) >=60) && (L1_f.sza(i_Lu) <70);
#          ksun0 = 2.47; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 2.03; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       elseif (L1_f.sza(i_Lu) >=70);
#          ksun0 = 2.33; # Point sensor from Gordon and Ding 1992 as per Mueller 2003 
#          ksun1 = 1.91; # Finite sensor value from Gordon and Ding 1992 as per Mueller 2003 
#       endif
# 
#       ksun(i_Lu) = ((1-g)*ksun0)+(g*ksun1);
#  
#    # Calculate K'sky
#       ksky(i_Lu) = 4.61 - 0.87*g;
# 
#    # Calculate E
#       for i_wv = 1:size(L2.instr.Lu.data,2);
#          esun(i_Lu,i_wv) = 1 - exp(-ksun(i_Lu)*int_a(i_Lu,i_wv)*0.25);
#          esky(i_Lu,i_wv) = 1 - exp(-ksky(i_Lu)*int_a(i_Lu,i_wv)*0.25);
#          etot(i_Lu,i_wv) = (esun(i_Lu,i_wv) + esky(i_Lu, i_wv))./(1+sf_h(i_Lu, i_wv));
#       endfor
# 
#    # Apply correction to Lu data
#       L2.Lucorr.data(i_Lu,:) = L2.instr.Lu.data(i_Lu,:)./(1-etot(i_Lu,:));
#    endfor
# 
#    ## sea-surface transmittance coefficient; Tzortziou et al., 2006 ##
#    L2.Lwcorr.data = L2.Lucorr.data * 0.544;
#    L2.Lw.data = L2.instr.Lu.data * 0.544;
# else
# all the above is commented out using the syntax "#{" and "#}"

 disp("----- skipped commented block -----");
 fflush(stdout);

% glint correction to water leaving radiance equation
   L2.rho_fitted = L2.rho; % rho_fitted is name for opean ocean rho method
   if size(L2.rho_fitted,2)==1; % 1 column implies aeronet rho method has been used
      L2.Lw.data = L2.instr.Lt.data - L2.rho_fitted.*L2.instr.Li.data;
   else % 2 columns implies open ocean rho method has been used
      L2.Lw.data = L2.instr.Lt.data - L2.rho_fitted(:,1).*L2.instr.Li.data - L2.rho_fitted(:,2);
   endif
   
#endif ??

### Apply filter for negative data ###
ikeep = find(L2.Lw.data(:,L2.wv==600)>0);
L2 = hsas_filter(L2, ikeep);

if strcmp(INSTRUMENT,'triosIW');
   ikeep = find(L2.Lw.data(:,L2.wv==600)<1);
   L2 = hsas_filter(L2, ikeep);
   ikeep = find(L2.Lw.data(:,L2.wv==400)>1);
   L2 = hsas_filter(L2, ikeep);
endif

### Compute Rrs ###
if strcmp(INSTRUMENT,'triosIW');
   L2.Rrs.data = L2.Lw.data./L2.instr.Es.data;
   L2.Rrs_corr.data = L2.Lwcorr.data./L2.instr.Es.data;

else
	L2.Rrs.data = L2.Lw.data./L2.instr.Es.data;

endif

% BDRF correction 
L2 = hsas_R0_R(L2); % 'gothic R' ratio for flat surface - see Gordon 2005 (assumes dependence on theta_v but not windspeed)


# Iteratively estimate chl or add ACS Chl to L2 structure and retrieve fQ values or use L2.chl to extract fQ
L2.conf.chl_alg = "A";  # (A = Atlantic)
L2.conf.used_acs = false; # flag: ACS Chl used


L2 = hsas_extract_fQ(L2, L2.conf.chl_alg);


# Calculate normalised radiometric quantities
L2 = hsas_extract_F0(L2);

keyboard
# Calculate Lw0
#L2.Lw.data,L2.R0_R,L2.fQs0_fQs
L2.Lw0.data = L2.Lw.data .* L2.R0_R .* L2.fQs0_fQs ;

# Calculate normalised water leaving radiance 
L2.exLwn.data = L2.Lw.data .* L2.R0_R .* L2.fQs0_fQs .* L2.fQ0_fQs0 .*   ((ones(size(L2.gps.lat)))*L2.F0 ./ L2.instr.Es.data);
L2.exRrs.data = L2.Lw.data .* L2.R0_R .* L2.fQs0_fQs .* L2.fQ0_fQs0 .* 1 ./ L2.instr.Es.data;

# Save file
# disp(["\n\n writing output file ", [fnout FILT_TAG ".mat \n\n"]]);
# fflush(stdout);

L2.files.L2 = fnout ;

save("-v7", fnout, "L2", "L1_f");
disp(["Written file " fnout "\n\n\n\n"]);
fflush(stdout);







% # write NetCDF4 file
[ tmp0, tmpfn] = fileparts(L2.files.L2);
fn_L2_NetCDF = [dout_nc tmpfn ".nc"];
hsas_write_L2_NetCDF(L2, fn_L2_NetCDF);








