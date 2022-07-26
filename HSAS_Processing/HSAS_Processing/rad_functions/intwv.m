function intwv = hsas_intwv(sn, in, wv, VERBOSE)


# interpolate to common wavelength range
    if VERBOSE, disp(["interpolating wavelengths of " sn{2} " spectra..."]); fflush(stdout); endif



# interpolate wavelengths

    intwv.data = interp1(in.wv, in.data', wv)';

    intwv.time = in.time;
    


    if isfield(in, 'int_time_sec')
        intwv.int_time_sec = in.int_time_sec;
    endif


    
    if isfield(in, 'int_time_sec')
        intwv.cal_file = in.cal_file;
    endif




endfunction
