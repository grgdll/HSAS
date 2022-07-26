function L2 = rd_true_wind(L2)

global DIR_SHIP    

   # read wind speed from anemometer
        fn_ws = [DIR_SHIP "anemometer/anemometer." L2.doy];
        ws = rd_anemometer(fn_ws);  

      # check if there are large gaps in the data (in which case the interpolation below should be replaced with an intersection
           max_delta_seconds_ws = max(diff(ws.time)*1440*60); # this is the max time difference between measurements
         if max_delta_seconds_ws > 60
               disp(["WARNING: there is a LARGE GAP in the data from the anemometer." L2.doy " file: replace interpolation with intersection (or something else)"] ); 
               fflush(stdout);
               L2.WARNINGS.(['anemometer_' L2.doy]) = ["WARNING: there is a LARGE GAP in the data from the anemometer." L2.doy " file: replace interpolation with intersection (or something else)"];
         endif


    # read vessel course and speed over ground
        fn_sog = [DIR_SHIP "seatex-vtg/seatex-vtg." L2.doy];
        vs = rd_seatex_vtg(fn_sog);

      # check if there are large gaps in the data (in which case the interpolation below should be replaced with an intersection
           max_delta_seconds_vs = max(diff(vs.time)*1440*60); # this is the max time difference between measurements
         if max_delta_seconds_vs > 60
               disp(["WARNING: there is a LARGE GAP in the data from the seatex-vtg." L2.doy " file: replace interpolation with intersection (or something else)"] ); 
               fflush(stdout);
               L2.WARNINGS.(['seatex_vtg_' L2.doy]) = ["WARNING: there is a LARGE GAP in the data from the seatex-vtg." L2.doy " file: replace interpolation with intersection (or something else)"];
 
         endif



    # interpolate wind data to radiometric data
        flds = fieldnames(vs);
        for ifld = 1:length(flds)
            if strcmp(flds{ifld}, "time")
                continue
            endif
            vs.ii.(flds{ifld}) = interp1(vs.time, vs.(flds{ifld}), L2.time);
        endfor
        
        flds = fieldnames(ws);
        for ifld = 1:length(flds)
            if strcmp(flds{ifld}, "time")
                continue
            endif
            ws.ii.(flds{ifld}) = interp1(ws.time, ws.(flds{ifld}), L2.time);
        endfor

    # compute true wind
    # jcr_truew(num, sel=3, vs.course_degs, tmp.sog_m2s, ws.wdir, zlr=0, L2.gps.hdg, adir, ws.ws, wmis=[], tdir, tspd, nw, nwam, nwpm, nwf)
       
    [apparent_wind_dir, true_wind_dir, true_wind_spd, nw, nwam, nwpm, nwf] = jcr_truew(num=length(vs.ii.course_degs), sel=3, vs.ii.course_degs, vs.ii.sog_m2s, ws.ii.wdir, zlr=0, L2.gps.hdg, ws.ii.ws, wmis=nan(num,5) );

    L2.wdir = true_wind_dir(:);
    L2.ws = true_wind_spd(:);    
    L2.sog_m2s = vs.ii.sog_m2s(:);  # [m/s]  


endfunction

