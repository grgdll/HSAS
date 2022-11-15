# read hsas digital counts and apply calibration coefficients


% Modifications for FICE 2022 by tjor:

# (i)  Calibration applies to stations rather than days.
# (ii) Non-linearity correction now implemented in correct_non_linearity_at_Cal_FICE.m (accomodates new data format from Tartu).
# (iii) Straylight correction has been modified to take Linear Spread Function matrix (L), from which straylight matrix D_SL in derived
# previous input was via D_SL directly)
 
 
clear all
close all

 
% pkg load financial

graphics_toolkit("gnuplot");

warning off #Turn off warnings 


addpath(strcat(pwd, "/cruise_specific_functions")) %assumes code is run from ../HSAS_Processing
addpath(strcat(pwd, "/rad_functions/"))
addpath(strcat(pwd, "/rad_functions/intwv"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa"))

# read input parameters for this cruise
input_parameters_hsas;

## these are the dirs containig the PRE and POST cals 
	if isempty(DIN_CALS_POST) 
	    din_cals = {DIN_CALS_PRE}; 
	else
	    din_cals = {DIN_CALS_PRE, DIN_CALS_POST}; 
	endif

## compare pre- and post-cruise cals
for iSN = 1:length(sn)
	
	disp([ "processing SN = " num2str(iSN) ]);
	fflush(stdout);
     
	 
	#initialize variables
	clear fncal cal offset_ gain_ wv_ int_time_
	 
    for ical = 1:length(din_cals)

        % fncal = glob([DIR_CAL din_cals{ical} "*" SN{iSN} "*.cal"]);
        fncal = glob([din_cals{ical} "*" sn{iSN} "*.cal"]);
        fncal_nl = glob([DIN_Non_Linearity NL_files_pre{ical}]); # fncal_nl 


        if length(fncal) > 2 # this is for when there are more cal files for one instrument (e.g., CAL_G SN223)
            fncal = sort(fncal){1};
        elseif length(fncal)==1   # tjor - else statement modified for FICE - I think we just want to read pre/post in turn for each 	 sensor?
            fncal = fncal{1};
        endif
        

        # read and store cal files
        cal{ical} = hsas_rd_satlantic_cal(fncal); # NL: correction 
        offset_(ical,:) = cal{ical}.offset;
        gain_(ical,:) = cal{ical}.gain;
        wv_(ical,:) = cal{ical}.wv;
        int_time_(ical,:) = cal{ical}.int_time_wv;       

    endfor
          
    # plot  post/pre-1 cal coefficients    
	if length(cal)==2    

       figure(1, 'visible', 'off');
        clf
        hold on
        subplot(121)
            plot(cal{1}.wv, cal{2}.offset./cal{1}.offset-1, [";" datestr(cal{1}.date, "yyyy/mmm/dd") "\n" datestr(cal{2}.date, "yyyy/mmm/dd") ";"])
            % ylim([-1 1]*mean(abs(cal{2}.offset./cal{1}.offset-1))*1.5) % commented due to nan padding
            set(gca, 'ygrid', 'on', 'gridlinestyle', ':');
            hold on, plot(cal{1}.wv, cal{1}.wv*0, 'k')
            xlim([350 850])
            
            xlabel('wavelength [nm]')
            ylabel('POST/PRE -1 ')
            title('offset')
            
        subplot(122)
            plot(cal{1}.wv, cal{2}.gain./cal{1}.gain-1, [";" datestr(cal{1}.date, "yyyy/mmm/dd") "\n" datestr(cal{2}.date, "yyyy/mmm/dd") ";"])
            ylim([-1 1]*0.05)
            set(gca, 'ygrid', 'on', 'gridlinestyle', ':');
            hold on, plot(cal{1}.wv, cal{1}.wv*0, 'k')
            xlim([350 850])
               
            xlabel('wavelength [nm]')
            #ylabel('ppost/pre -1 [-]')
            title('gain')

        
        set(gcf, 'paperposition', [0.25 0.25 12 4])
        fnout = [DIR_CAL "/rel_change_in_cal_coeffs_" strsplit(fncal, {"/","."}){end-1}(1:end-1) ".png"]   ;
        print("-dpng", fnout);
	endif 
    
      
    ## compute and write average cal file (+ uncertainty) for this instrument SN
        # compute stats of cal coeffs
		if length(cal)==2
            gain.mean = mean(gain_);
            gain.unc = std(gain_);
            gain.n = size(gain_,1);
            
            offset.mean = mean(offset_) ;
            offset.unc = std(offset_);
            offset.n = size(offset_,1);
            
            int_time = mean(int_time_);
            
            if ical > 1 & (int_time_(1,20) != int_time_(2,20)) % ~all(std(int_time_,[],1)<=eps)
                disp('integration time has changed between calibrations!!!');
                keyboard
            endif
            
            wv = mean(wv_);
            
            
            if ical>1 & (wv_(1,20) != wv_(2,20))%~all(std(wv_,[],1)<=eps)
                disp('wavelengths have changed between calibrations!!!');
                keyboard
            endif
			
		elseif length(cal)==1
			
            gain.mean = gain_;
            gain.unc = gain_*nan;
            gain.n = size(gain_,1);
            
            offset.mean = offset_ ;
            offset.unc = offset_*nan;
            offset.n = size(offset_,1);
            
            int_time = mean(int_time_);
  
            wv = wv_;
	
		endif	
				
			
	# save existing calibrations in octave binary format
       	if ~exist([DIR_CAL CRUISE])    
            mkdir([DIR_CAL CRUISE]);  
        endif
        
       	fnout_cal = [DIR_CAL CRUISE "/mean_" radiometers{iSN} sn{iSN} ".cal"]   ;    
       	save("-mat-binary", fnout_cal, "cal", "gain", "offset", "wv", "int_time");
        
   
	 	
 % Read serial number and sensor type
	rad_sn = cell2struct(sn,radiometers,2);
        sn_rad = cell2struct(radiometers, sn, 2);
        sensor_id = sn{iSN};	

# tjor: modifications for FICE non-linearity follow from this point :	
# FICE deployment has new data format for NL coefficients; see `radcal_recipies.txt' in HyperSAS_config for more details
# We now read a general data matrix, from which alpha is calculated within `correct_non_linearity_at_Cal_FICE()'
# We also average pre and post deployment data martrices (similar approach to cal files), before deriving alpha.
 if FLAG_NON_LINEARITY == 1
	  pkg load io
	  #---radiometer related to sn
	  	  disp('Loading Non-linearity correction coefficients....')  
	 
	  # pre deployment
        data_range = [1345,0,1525,10]; # reads new data matrix for each sensor callibration
	coeff_ES_pre = dlmread([DIN_Non_Linearity NL_files_pre{1}],'', data_range);
        data_range = [1484,0,1664,10];
	coeff_LI_pre = dlmread([DIN_Non_Linearity NL_files_pre{2}],'', data_range);
	data_range = [1484,0,1664,10];
	coeff_LT_pre = dlmread([DIN_Non_Linearity NL_files_pre{3}],'', data_range);
	  
	  # post deployment
	data_range = [115,0,295,10]; # reads new data matrix for each sensor callibration
	coeff_ES_post = dlmread([DIN_Non_Linearity NL_files_post{1}],'', data_range);
	data_range = [254,0,434,10];
	coeff_LI_post = dlmread([DIN_Non_Linearity NL_files_post{2}],'', data_range);
	data_range = [254,0,434,10];
	coeff_LT_post = dlmread([DIN_Non_Linearity NL_files_post{3}],'', data_range);

	# plot function for 
        # figure(1, 'visible', 'off');
        #  clf
        #  subplot(131)
        #  plot(coeff_ES_pre(:,1), coeff_ES_pre(:,2))
        #  hold on
        #  plot(coeff_ES_post(:,1), coeff_ES_post(:,2),'r--')
        #  legend ("Pre", "Post");
        #  % ylim([-1 1]*mean(abs(cal{2}.offset./cal{1}.offset-1))*1.5) % commented due to nan padding
        #  set(gca, 'ygrid', 'on', 'gridlinestyle', ':');
        #  xlim([350 850])
        #  ylim([0,0.001])
        #  xlabel('wavelength [nm]')
        #  ylabel('cal coefficient ')
        #  title('Es: Non-linearity cal coefficients (pre / post')
	  
	#  subplot(132)
        #  plot(coeff_LI_pre(:,1), coeff_LI_pre(:,2))
        #  hold on
        #  plot(coeff_LI_post(:,1), coeff_LI_post(:,2),'r--')
        #  legend ("Pre", "Post");
        #  % ylim([-1 1]*mean(abs(cal{2}.offset./cal{1}.offset-1))*1.5) % commented due to nan padding
        #  set(gca, 'ygrid', 'on', 'gridlinestyle', ':');
        #  xlim([350 850])
        #  ylim([0,0.001])
        #  xlabel('wavelength [nm]')
        #  ylabel('cal coefficient ')
        #  title('Li: Non-linearity cal coefficients (pre / post')
	
	#  subplot(133)
        #  plot(coeff_LT_pre(:,1), coeff_LT_pre(:,2))
        #  hold on
        #  plot(coeff_LT_post(:,1), coeff_LT_post(:,2),'r--')
        #  legend ("Pre", "Post");
        #  % ylim([-1 1]*mean(abs(cal{2}.offset./cal{1}.offset-1))*1.5) % commented due to nan padding
        #  set(gca, 'ygrid', 'on', 'gridlinestyle', ':');
        #  xlim([350 850])
        #  ylim([0,0.001])
        #  xlabel('wavelength [nm]')
        #  ylabel('cal coefficient ')
        #  title('Li: Non-linearity cal coefficients (pre / post)')
  
	 
	#  set(gcf)
	#  fnout = [DIR_CAL "/NLcoeff_pre_post.png"]   ;
	#  print("-dpng", fnout);
	 
	 
	# average of pre/post data matrices for each sensor callibration
           coeff_ES =  0.5*(coeff_ES_pre + coeff_ES_post);
           coeff_LI =  0.5*(coeff_LI_pre + coeff_LI_post);
	   coeff_LT =  0.5*(coeff_LT_pre + coeff_LT_post);
	 
	   non_linearity_coeff = struct('coeff_LI',coeff_LI,'coeff_LT',coeff_LT,'coeff_ES',coeff_ES);
  else
  	   non_linearity_coeff = struct('coeff_LI',nan,'coeff_LT',nan,'coeff_ES',nan);
  endif
    
		
  #----Read Straylight Distribution Matrix
  # tjor - for 2022 SL files, Tartu now provide the linear spread function matrix (LSF)
  # This is now normalized to give the straylight distribution matrix, D_SL. eg. see,  
  # https://iopscience.iop.org/article/10.1088/0952-4746/34/4/915
  # 
   if FLAG_STRAY_LIGHT == 1
	  disp('Loading StrayLight Distribution Matrix....')  

	  if sn_rad.(sensor_id) == 'ES'
	 	fn = [DIR_SLCORR, FN_SLCORR_ES];
	  endif
	  if sn_rad.(sensor_id) == 'LI'
	  	fn = [DIR_SLCORR, FN_SLCORR_LI];
	  endif
	  if sn_rad.(sensor_id) == 'LT'
	  	fn = [DIR_SLCORR, FN_SLCORR_LT];
	  endif
	  
	  LSF = dlmread(fn,'',[49,0,304,255]); # read LSF (linear spread function) is 256 X 256 matrix, which matches wl bins 
	  # in rad_cal files 	 
	 
	  # convert to straylight matrix D_SL
	  D = LSF./(sum(LSF')); # out-of-band normalization  - normalize w.r.t. sum over each row of LSF matrix (energy of each scan)
	  for i=1:length(D)    
	  	D(i,i) = 0;     # in-band normalization - set diagonal elements to zero
	  endfor 
	
	  D_SL = D(cal{1}.usedpixels,cal{1}.usedpixels); # Reduces matrix to 137 * 137 (matching wl bin indices in rad_cal 		  		  files) - bin 14 is normally the first wl bin used in processing, 150 is the last #
	  
  else
  	  D = D_SL = nan;
 endif
  

		
	###### calibrate data from this instrument
	days = glob([DIN_HSAS "2022*"]); # read all days (stations for fice 2022)
	istart = min(find(cellfun(@isempty, strfind(days, DAY_START))==0)); # find first day (station) to be processed
	istop = max(find(cellfun(@isempty, strfind(days, DAY_STOP))==0)); # find last day (station) to be processed
	
	
	if isempty(istart) | isempty(istop)
		disp('cannot find start or stop day');
		keyboard
	endif

 
	
	# loop over the days that need to be processed
	for iday = istart:istop
		
		sday = strsplit(days{iday}, '/'){end}; # extract date string 
		# these are the files to process
        	fn = glob([DIN_HSAS sday "/*" sn{iSN} "*.dat"]);

	        for ifn = 1:length(fn)
        
				disp(sprintf("%u/%u", ifn, length(fn)));
				fflush(stdout);

		        # read digital counts
		
	            	out = hsas_rd_digital_counts(fn{ifn});
				    if isempty(out.time)
						continue;
				    endif            
			
	      # tjor - extra filter step for FICE. This is done so that `out' corresponds to same wl bins as offset. Ideally, in future, 	processing, this should be done via modifying hsas_rd_digital counts. 
	           if length(out.wv) > length(cal{1}.wv)
	           	out.(out.instru) = out.(out.instru)(:,cal{1}.usedpixels)
	           	out.wv = cal{1}.wv
	           endif
	  
	                ff = [out.instru "cal"];
	                out.cal_file = fnout_cal;
                    out.(ff) = hsas_calibrate_with_correction(sn{iSN}, out.wv, L_CountsLightDat=out.(out.instru), L_CalDarkDat=offset.mean, gain.mean,immers_coeff=1, it_1=int_time, it_2=out.int_time_sec, rad_sn, sn_rad, non_linearity_coeff,D_SL,FLAG_NON_LINEARITY,FLAG_STRAY_LIGHT);
	                #out.(ff) = hsas_calibrate(L_CountsLightDat=out.(out.instru), L_CalDarkDat=offset.mean, gain.mean, immers_coeff=1, it_1=int_time, it_2=out.int_time_sec);

	            ### sort fields alphabetically
	                out = orderfields(out);
                
	            ### write calibrated files (which format?)       
					dout = [DOUT_HSAS "Calibrated/" sday "/"]; # create dir for calibrated files
					if ~exist(dout)
						mkdir(dout);
					endif
					

	                fnout = strrep(out.satcon_file, "Raw/RawExtracted", [DOUT_HSAS_SUB "Calibrated"]);
	                hsas_write_calibrated_Satlantic_file_format(fnout, out, DELIM=' ');
                	
					disp(["Written calibrated file: " fnout]);
					fflush(stdout);
					
	        endfor # fn
      
	  
	endfor # days
	    

endfor # sn













