function cal = hsas_calibrate_with_correction(sensor_id, wl, L_CountsLightDat, L_CalDarkDat, a, ic, it_1, it_2, rad_sn, sn_rad, non_linearity_coeff, D_SL,FLAG_NON_LINEARITY,FLAG_STRAY_LIGHT)
% cal = hsas_calibrate(L_CountsLightDat, L_CalDarkDat, a, ic, it_1, it_2)
#
# cal = (L_CountsLightDat - L_CalDarkDat).*a.*ic.*it_1./it_2
#	
# where:
# # a is the calibration coefficient, 
# # ic is an immersion coefficient, 
# # it_1 is the integration time during calibration and 
# # it_2 is the integration time during the measurement. 
#
# a, ic and it1 are taken from a calibration file, it2 is obtained from the same log file as optical data.
	
	#---------------------------------



################################# Correction ######################


data = (L_CountsLightDat - L_CalDarkDat);




#---- nonlinear correction - #
if FLAG_NON_LINEARITY == 1
	data = correct_non_linearity_at_Cal_FICE(rad_sn, sn_rad, sensor_id, non_linearity_coeff,data,wl);
endif



#---Stray Light Correction
if FLAG_STRAY_LIGHT == 1
	data = hsas_straylight_correct_at_Cal(sensor_id, wl, sn_rad, D_SL, data);
endif


cal = data.*a.*ic.*it_1./it_2;

	
	
endfunction
