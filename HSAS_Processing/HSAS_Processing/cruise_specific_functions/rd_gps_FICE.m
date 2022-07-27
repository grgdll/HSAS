function gps = rd_gps_FICE(fn)
% This is a modified GCPS function for the FICE 2022 deployment - it will require replacement of read_seatex_gga

# function gps = rd_JCR_gps(fn)
# read gps data from JCR's underway data stream
# this is thought to be similar to rd_DY_gps.m, but neews to read two different files for JCR
#
# fn is the path to the seatex-gga.ACO file

	input_parameters_hsas;

	pkg load io
	keyboard
	DATESTR = strsplit(fn, {'/','_'}){12};
    	[Data, Headers] = xlsread(fn)

	gps = rd_seatex_gga(  fn   );

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

	
	
