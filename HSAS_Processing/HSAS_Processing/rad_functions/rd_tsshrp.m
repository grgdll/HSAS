function tmp = rd_tsshrp(fn)

#2016,297.000000,297,0.00000000,7.694118    ,-1043.854887,-4.000000,U     ,-172.000000,144.000000
#yyyy,doy.decdoy,doy,decdoy    ,hacc[cm/s^2],vacc[cm/s^2],heave[cm],status,roll[degs] ,pitch[degs]

#### first make sure there are no double commas in the file
#    system(["rm -f tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
#    system(["sed s/,,/,nan,/ " fn " > tmp" strsplit(fn ,"/"){end-1} ".tmp"]);
#    system(["mv tmp" strsplit(fn ,"/"){end-1} ".tmp " fn]);





    fmt = "%f,%f,%f,%f,%f,%f,%f,%*c,%f,%f\n";

    fid = fopen(fn, "r");

        d = fscanf(fid, fmt, [9, Inf])'; 
   
    fclose(fid);   


    y0 = datenum(d(1),1,1);
    tmp.time = y0-1+d(:,2);
    tmp.roll = d(:,end-1)/100; # [degrees]
    tmp.pitch = d(:,end)/100; # [degrees]
 
    tmp.tilt = hsas_cmp_tilt(tmp.roll, tmp.pitch);




endfunction
