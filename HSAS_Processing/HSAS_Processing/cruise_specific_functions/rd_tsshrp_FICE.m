function tmp = rd_tsshrp_FICE(DATESTR)
% version of ATT function for FICE 2022 (no tilt/roll data)
	


    tmp.time = datenum(DATESTR,'yyyymmdd') + [0.2, 0.8]; % fake vector for time interpolaton
    tmp.heave = [0 ,0]; # [cm]
    tmp.roll = [0, 0]; # [degrees]
    tmp.pitch = [0, 0]; # [degrees]
 

endfunction
