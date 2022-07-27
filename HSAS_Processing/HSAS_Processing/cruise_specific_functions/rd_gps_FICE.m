function gps = rd_gps_FICE(fn)
% This is a modified GPS function for the FICE 2022 deployment - it reads meta data from an xls spreadsheet.
% Lat-lon are constant; data fields not present (but are usually in ship GPS files) are assigned NaN
% The relative solar azimuth (phi) is read directly


	input_parameters_hsas;

	pkg load io

	DATESTR = strsplit(fn, {'/','_'}){12};
    	[data, headers] = xlsread(fn, DATESTR); % DATESTR is used to pass sheet number
    	station_number = data(:,1); 
    	id = headers(2:end,2); % Identifier for metadata record
    	orientation = headers(2:end,3);
    	phi = nan(length(orientation),1); % phi is relative (signed) sensor azimuth w.r.t. sun
    	for ist = 1:length(orientation)
    		tmp = strsplit(orientation{ist},' ');
    		phi(ist) = str2num(tmp{1});
    		if strcmp(tmp{2}, 'sx') % sensor with azmiuth to left of sun defined as -ve, right azimuth as +ve
    			phi(ist) = -phi(ist);
    		endif 
    	endfor
    		
	keyboard
	# gps = rd_seatex_gga(  fn   );

	# now read speed and course over ground
	fn2 = strrep(fn, "gga", "vtg");
	tmp = rd_seatex_vtg(  fn2   );
	
	# interpolate tmp onto gps
	gps.sog_m2s = interp1(tmp.time, tmp.sog_m2s, gps.time);
	gps.cog_deg = interp1(tmp.time, tmp.course_degs, gps.time);
	

	# now read heading
	fn3 = strrep(fn, "gga", "hdt");
	tmp2 = rd_seatex_hdt(  fn3   );
	
	# interpolate tmp2 onto gps
	gps.hdg = interp1(tmp2.time, tmp2.hdg, gps.time);
	
	
endfunction

	
	
