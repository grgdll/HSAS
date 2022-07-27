function gps = rd_gps_FICE(fn)
% This is a modified GPS function for the FICE 2022 deployment - it reads meta data from an xls spreadsheet.
% Lat-lon are constant; data fields not present (but are usually in ship GPS files) are assigned NaN
% The relative solar azimuth (phi) is read directly
% This replaces output of rd_seatex


	input_parameters_hsas;

	pkg load io

	DATESTR = strsplit(fn, {'/','_'}){12}; 
    	[data, headers] = xlsread(fn, DATESTR); % DATESTR is used to pass sheet number
    	station_number = data(:,1); 
    	id = cell2mat(headers(2:end,2)); % Identifier for metadata record
    	orientation = headers(2:end,3);
    	
    	L = length(orientation);
    	phi = nan(L,1); % phi is relative (signed) sensor azimuth w.r.t. sun
    	for ist = 1:L
    		tmp = strsplit(orientation{ist},' ');
    		phi(ist) = str2num(tmp{1});
    		if strcmp(tmp{2}, 'sx') % sensor with azmiuth to left of sun defined as -ve, right azimuth as +ve
    			phi(ist) = -phi(ist);
    		endif 
    	endfor
    		
    	% fill gps structure
	gps.time = datenum(id(:,21:end),'yyyymmdd_HHMMSS'); 
	gps.lat = 45.31425*ones(L,1); % FICE coordinates (hardcoded step)
	gps.lon = 12.508317*ones(L,1);
  	gps.phi = phi;
  	gps.sog_m2s = nan(L,1);
  	gps.cog_deg = nan(L,1);
  	gps.hdg = nan(L,1);

 
	
endfunction

	
	
