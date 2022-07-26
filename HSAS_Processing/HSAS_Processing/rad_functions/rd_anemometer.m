function d = rd_anemometer(fn)


addpath ~/Dropbox/Octave_functions

#### first make sure there are no double commas in the file
    system(["rm -f tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
    system(["sed s/,,/,nan,/ " fn " > tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
    system(["mv tmp" strsplit(fn ,"/"){end-1} ".tmp " fn]);

    fid = fopen(fn, "r");

#  amt24         2014,268.999410,268,0.99940972,005,18.5
        fmt = "%f,%f,%f,%f,%f,%f,%f\n";
        nFields = length(find(cellfun(@isempty, strfind(strsplit(fmt, ","), "\*"))==1));

        tmp = fscanf(fid, fmt, [nFields, inf])';

    fclose(fid);

    y0 = datenum(tmp(1),1,1);

    d.time = y0-1+tmp(:,2);

    d.ws = tmp(:,end); # [m/s]
    
    d.wdir = tmp(:,end-1) ; # [degrees]


endfunction



