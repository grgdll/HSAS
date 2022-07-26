function L2 = hsas_filter_sensors_using_Lu_data_v2(L1_f, L2, vza_nm)


# this is because I changed the name of the variable in the structure that contains the vza
    if length(argn)==3
        vza_nm = "vza";
    endif



##### filter : blo changed 850 to 820 to avoid nan's at the hig NIR end wavelengths
    [tmp iNIR] = intersect(L2.wv, 700:820); # find index for NIR wavelengths 
    x = mean(L1_f.instr.Lu.data(:,iNIR)');; # take the mean of the NIR data 

    delta_t_secs = round((median(diff(L1_f.time)))*60*60*24);  # [seconds]  This is the median time-difference between scans
    if isfield(L2.instr.Es, "sn")   # this is for Trios
#        time_window = 15;  # data for trios are every 4 seconds so this corresponds to a 1-minute window
        time_window = 60/delta_t_secs;  # computes approx how many scans we have in one minute
    else    # this is for Hsas
        time_window = 60/delta_t_secs;  # computes approx how many scans we have in one minute
    endif
    time_window = round(time_window);
    R = slidefun('min', time_window, x,'central');
    iLu = find(x==R)';
    L2.iLu = iLu;



    
    
##### apply filter to Lu data 
    flds = fieldnames(L1_f.instr);

    for i_instr = 1:length(flds);

        # CHANGE HERE WHEN WE WANT TO FILTER ALL SENSORS   
#        if strcmp(flds{i_instr}, "Lu")      # COMMENT THIS
            ifilt = iLu;                     # DO NOT COMMENT THIS
#        else                                # COMMENT THIS
#            ifilt = 1:length(L1_f.lat);;    # COMMENT THIS
#        endif                               # COMMENT THIS

        L2.instr.(flds{i_instr}).data = L1_f.instr.(flds{i_instr}).data(ifilt,:);

        if isfield(L1_f.instr.(flds{i_instr}), 'pitch');
            L2.instr.(flds{i_instr}).pitch = L1_f.instr.(flds{i_instr}).pitch(ifilt,:);
            L2.instr.(flds{i_instr}).roll = L1_f.instr.(flds{i_instr}).roll(ifilt,:);
            L2.instr.(flds{i_instr}).tilt = L1_f.instr.(flds{i_instr}).tilt(ifilt,:);
            L2.instr.(flds{i_instr}).(vza_nm) = L1_f.instr.(flds{i_instr}).(vza_nm)(ifilt,:);

        else
            L2.instr.(flds{i_instr}).pitch = L1_f.pitch(ifilt,:);
            L2.instr.(flds{i_instr}).roll = L1_f.roll(ifilt,:);
            L2.instr.(flds{i_instr}).tilt = L1_f.tilt(ifilt,:);
            L2.instr.(flds{i_instr}).(vza_nm) = L1_f.(vza_nm)(ifilt,:);


        endif


        if isfield(L1_f.instr.(flds{i_instr}), 'int_time_sec')
            L2.instr.(flds{i_instr}).int_time_sec = L1_f.instr.(flds{i_instr}).int_time_sec(ifilt,:);
# blo: commented out as this variable is NOT in trios, so we don't want to add it if performing two filters
#        else 
#            L2.instr.(flds{i_instr}).int_time_sec = "";
        endif

        if isfield(L1_f.instr.(flds{i_instr}), 'time')
            L2.instr.(flds{i_instr}).time = L1_f.instr.(flds{i_instr}).time(ifilt,:);
        else 
            L2.instr.(flds{i_instr}).time = L1_f.time(ifilt,:);
        endif
#        L2.instr.(flds{i_instr}) = hsas_cmp_stats(L2.instr.(flds{i_instr}));

    endfor



    flds = fieldnames(L2);

#    filter all data except for instrument data   
        for i_fld = 1:length(flds);

            if isstruct(L2.(flds{i_fld})) || size(L2.(flds{i_fld}),1)!=size(L1_f.gps.lat,1);
                continue
            endif

            L2.(flds{i_fld}) = L2.(flds{i_fld})(iLu,:);

        endfor

        
        L2.gps.time = L1_f.gps.time(ifilt,:);
        L2.gps.lat = L1_f.gps.lat(ifilt,:);
        L2.gps.lon = L1_f.gps.lon(ifilt,:);
        L2.gps.hdg = L1_f.gps.hdg(ifilt,:);


        L2.time = L2.instr.Es.time;


endfunction
