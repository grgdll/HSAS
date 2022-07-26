function tmp = rd_gps_discovery(in_folder, gps_date)
   % Read GPS variables from netcdf GPS file
   % date is used to identify the file

   pkg load octcdf
   
   ncfiles = glob([in_folder, gps_date '*position-Applanix_GPS*']);

   % There should be only one file returned by glob
   if length(ncfiles)~=1
      disp('Something wrong with GPS files')

   else
      ncfile = ncfiles{1};
   end
   nc =netcdf(ncfile,'r');
   tmp.time = nc{'time'}(:)+datenum([1899,12,30,0,0,0]);
   tmp.lon = nc{'long'}(:);
   tmp.lat = nc{'lat'}(:);
   tmp.hdg = nc{'heading'}(:);

   endfunction
