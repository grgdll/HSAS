function out = hsas_int2Lt(in, Lt_time)
# interpolate to Lt time step


# interpolate wavelengths

    out.data = interp1(in.time, in.data, Lt_time);
    
#    if isfield(in, 'int_time_sec')
#        out.int_time_sec = interp1(in.time, in.int_time_sec, Lt_time);

    out.time = Lt_time;
    
    if isfield(in , 'cal_file')
        out.cal_file = in.cal_file;
    else 
        out.cal_file = "";
    endif

endfunction
