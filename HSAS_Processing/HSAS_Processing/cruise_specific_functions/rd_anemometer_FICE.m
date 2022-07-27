function d = rd_anemometer_FICE(fn)
% FICE version of function to read wind data (xlsx files)
% Other metadata fields are available in xls meteo file.

    pkg load io
    [data, headers] = xlsread(fn);

    d.time = data(:,1) + 693960 - 1/24; % [decimal days] 693960 converts between excel and octave time, -1/24 accounts for UTC +1 to UTC

    d.wspd = data(:,4); # [m/s] - median in previous 10 mins
    
    d.wdir = data(:,3) ; # azimuth from N [degrees] - median in previous 10 mins


endfunction



