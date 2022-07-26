function d = rd_anemometer_FICE(fn)
%FICE version of function to read wind data (xlsx files)

    pkg load io
    [Data, Headers] = xlsread(fn)
    
    % note - now need to convert into structure d with fields time, wind speed and wind direction
    

#  amt24         2014,268.999410,268,0.99940972,005,18.5    
        fmt = "%f,%f,%f,%f,%f,%f\n";
        nFields = length(find(cellfun(@isempty, strfind(strsplit(fmt, ","), "\*"))==1));
        tmp = fscanf(fid, fmt, [nFields, inf])';

    fclose(fid);

    y0 = datenum(tmp(1),1,1);

    d.time = y0-1+tmp(:,2);

    d.wspd = tmp(:,end); # [m/s]
    
    d.wdir = tmp(:,end-1) ; # [degrees]


endfunction



