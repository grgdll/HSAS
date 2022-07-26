# read hasa digital counts and apply calibration coefficients
clear all
close all
 
% pkg load financial

graphics_toolkit("gnuplot");

warning off #Turn off warnings

addpath("../JCR")
addpath("../JCR/cruise_specific_functions")
addpath("../JCR/rad_functions/")
addpath("../JCR/rad_functions/intwv")
addpath("../JCR/rad_functions/DISTRIB_fQ_with_Raman")
addpath("../JCR/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa")


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
        
        if length(fncal)>2 # this is for when there are more cal files for one instrument (e.g., CAL_G SN223)
            fncal = sort(fncal){1};
        elseif
            fncal = fncal{2};
        endif
        
        # read and store cal files
        cal{ical} = hsas_rd_satlantic_cal(fncal);
        
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
            ylim([-1 1]*mean(abs(cal{2}.offset./cal{1}.offset-1))*1.5)
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
            
            if ical>1 & ~all(std(int_time_,[],1)<=eps)
                disp('integration time has changed beytween calibrations!!!');
                keyboard
            endif
            
            wv = mean(wv_);
            
            if ical>1 & ~all(std(wv_,[],1)<=eps)
                disp('wavelengths have changed beytween calibrations!!!');
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
        
        
		

		
	
  ####### Read non-linearity correction coefficients ##########
  pkg load io
  #---radiometer related to sn
  disp('Loading Non-linearity correction coefficients....')  
  rad_sn = cell2struct(sn,radiometers,2);
  sn_rad = cell2struct(radiometers, sn, 2);
  coeff_LI = xlsread(DIN_Non_Linearity,rad_sn.LI);
  coeff_LT = xlsread(DIN_Non_Linearity,rad_sn.LT);
  coeff_ES = xlsread(DIN_Non_Linearity,rad_sn.ES);
  non_linearity_coeff = struct('coeff_LI',coeff_LI(:,1:2),'coeff_LT',coeff_LT(:,1:2),'coeff_ES',coeff_ES(:,1:2));	
    
		
  #----Read Straylight Distribution Matrix
  disp('Loading StrayLight Distribution Matrix....')  
  sensor_id = sn{iSN};
  if sn_rad.(sensor_id) == 'ES'
  fn=[DIR_SLCORR,FN_SLCORR_ES];
  end
  if sn_rad.(sensor_id) == 'LI'
  fn=[DIR_SLCORR,FN_SLCORR_LI];
  end
  if sn_rad.(sensor_id) == 'LT'
  fn=[DIR_SLCORR,FN_SLCORR_LT];
  end

  D = load(fn);
  D_SL = D / norm(D);

		
	###### calibrate data from this instrument
    DIN_HSAS
	days = glob([DIN_HSAS "*"]); # read all days
	istart = find(cellfun(@isempty, strfind(days, DAY_START))==0); # find first day to be processed
	istop = find(cellfun(@isempty, strfind(days, DAY_STOP))==0); # find last day to be processed
	
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
			
	            ### apply cal
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













