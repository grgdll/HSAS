function tmp = rd_seatex_vtg(fn)

#2016,280.000000,280,0.00000000,165.55,,9.6,17.7


#### first make sure there are no double commas in the file
    system(["rm -f tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
    system(["sed s/,,/,nan,/ " fn " > tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
    system(["mv tmp" strsplit(fn ,"/"){end-1} ".tmp " fn]);





    fmt = "%f,%f,%f,%f,%f,%f,%f,%f\n";

    fid = fopen(fn, "r");

        d = fscanf(fid, fmt, [8, Inf])'; 
   
    fclose(fid);   

    y0 = datenum(d(1),1,1);

    tmp.time = y0-1+d(:,2);

    tmp.course_degs = d(:,end-3); # [degrees]
    tmp.sog_km2hr = d(:,end); # [km/hour]
    tmp.sog_m2s = tmp.sog_km2hr/3.6; # [m/s]



#    # this is to deal with change in format through the cruise
#    if length(tmp.time==1)
#
#        fmt = "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,,\n";
#
#        fid = fopen(fn, "r");
#
#            d = fscanf(fid, fmt, [11, Inf])'; 
#       
#        fclose(fid);   
#
#
#        tmp.time = y0(d(1))-1+d(:,2);
#
#        tmp.lat = d(:,6); # [degN]
#        tmp.lon = d(:,7); # [degE]
#            
#
#
#    endif




endfunction
