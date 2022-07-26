function cal = hsas_calibrate(L_CountsLightDat, L_CalDarkDat, a, ic, it_1, it_2)
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
	
	
	
	cal = (L_CountsLightDat - L_CalDarkDat).*a.*ic.*it_1./it_2;
	
	
	
endfunction
