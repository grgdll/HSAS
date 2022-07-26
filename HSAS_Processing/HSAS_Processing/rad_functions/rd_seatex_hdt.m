function tmp = rd_seatex_hdt(fn)

#2016,277.000012,277,0.00001157,194.21

    disp('reading seatex-hdt file...');
    fflush(stdout);


    fmt = "%f,%f,%f,%f,%f\n";

    fid = fopen(fn, "r");

        d = fscanf(fid, fmt, [5, Inf])'; 
   
    fclose(fid);   

    y0 = datenum(d(1),1,1);
    tmp.time = y0-1+d(:,2);

    tmp.hdg = d(:,end); # [degrees]
    





endfunction
