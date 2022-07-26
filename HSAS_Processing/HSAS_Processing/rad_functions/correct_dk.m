function nodk = hsas_correct_dk(sn, in, VERBOSE)
# interpolate and remove dark counts from raw data

    if VERBOSE, disp(["removing darks from " sn{2} " spectra..."]); fflush(stdout); endif

		
# interpolate dark readings
    isdk = find(in.dark);
    isnotdk = find(~in.dark);
    int_dk = interp1(in.time(isdk), in.data(isdk,:), in.time(isnotdk));


# correct for dark counts
    nodk.data = in.data(isnotdk,:)-int_dk;
    nodk.time = in.time(isnotdk);

    nodk.int_time_sec = in.int_time_sec(isnotdk);

    nodk.cal_file = in.cal_file;

    nodk.wv = in.wv;

    if in.sn == 223
        nodk.lat = in.lat(isnotdk);
        nodk.lon = in.lon(isnotdk);
    endif

endfunction





