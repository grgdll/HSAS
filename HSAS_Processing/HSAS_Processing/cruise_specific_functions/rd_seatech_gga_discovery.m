function tmp = rd_seatech_gga_discovery(date_str)
   % Read GPS variables from netcdf GPS file
   % date is used to identify the file

   global gps_dir

   ncfiles = glob([gps_dir, date_str '*position-Applanix_GPS*']);
   % There should be only one file returned by glob
   if length(ncfiles)~=1
      disp('Something wrong with GPS files')
      keyboard
   else
      ncfile = ncfiles{1};
   end%if
   nc =netcdf(ncfile,'r');
   % Time must be first element of tmp!!!
   % (otherwise error in step2h_underway_amt27_make_processed.m)
   % Assumes time is in dats (matlab format)
   tmp.time = nc{'time'}(:)+datenum([1899,12,30,0,0,0]);

   tmp.lon = nc{'long'}(:);
   tmp.lat = nc{'lat'}(:);

   %    % this is to deal with change in format through the cruise
   %    if length(tmp.time==1)
   %
   %        fmt = "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,,\n";
   %
   %        fid = fopen(fn, "r");
   %
   %            d = fscanf(fid, fmt, [11, Inf])'; 
   %       
   %        fclose(fid);   
   %
   %
   %        tmp.time = y0(d(1))-1+d(:,2);
   %
   %        tmp.lat = d(:,6); % [degN]
   %        tmp.lon = d(:,7); % [degE]
   %            
   %
   %
   %    end%if

   endfunction
