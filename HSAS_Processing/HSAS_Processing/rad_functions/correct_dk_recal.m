function nodk = hsas_correct_dk_recal(sn, in, VERBOSE)
# interpolate and remove dark counts from raw data

    if VERBOSE, disp(["removing darks from " sn{2} " spectra..."]); fflush(stdout); endif


# interpolate dark readings
    isdk = find(in.dark);
    isnotdk = find(~in.dark);
    int_dk = interp1(in.time(isdk), in.data(isdk,:), in.time(isnotdk));


# correct for dark counts
    nodk.data = in.data(isnotdk,:)-int_dk;

# apply recalibration scale factors after dark count has been removed
    nodk.data = nodk.data./in.recal_scaling(isnotdk,:);
    nodk.time = in.time(isnotdk);

    nodk.int_time_sec = in.int_time_sec(isnotdk);

    nodk.cal_file = in.cal_file;

    nodk.wv = in.wv;

    if in.sv == 223
        nodk.lat = in.lat(isnotdk);
        nodk.lon = in.lon(isnotdk);
    endif

endfunction





