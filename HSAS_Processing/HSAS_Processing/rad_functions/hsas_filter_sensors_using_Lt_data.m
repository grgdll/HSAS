function [L2] = hsas_filter_sensors_using_Lt_data(L1_f, L2, XX, vza_nm)


# this is because I changed the name of the variable in the structure that contains the vza
    if length(argn)==3
        vza_nm = "vza";
    endif



L2.threshold4Lt_filter = XX;
L2.conf.threshold4Lt_filter_in_percent = XX;


# no filter case
    if XX==100;

        iLt = 1:length(L1_f.instr.Lt.data(:,1));

        return
    endif

    







##### find lower XX% of Lt data (Marco fa cosi')

    [tmp iwv] = intersect(L2.wv, 450:650);

    tmp = mean(L1_f.instr.Lt.data(:,iwv)');

    iLt = find(tmp <= prctile(tmp',XX));
    
#    # metodo "casereccio" che usa Marco
#        [Stmp, Itmp] = sort(tmp(:), 1) ;
#        ii = floor(length(Itmp)/100*XX);
#        if ii==0, ii=1; endif    
#        iLt = Itmp(1:ii);


    
    
##### apply filter to Lt data 
    flds = fieldnames(L1_f.instr);

    for i_instr = 1:length(flds);

        # CHANGE HERE WHEN WE WANT TO FILTER ALL SENSORS   
#        if strcmp(flds{i_instr}, "Lt")      # COMMENT THIS
            ifilt = iLt;                     # DO NOT COMMENT THIS
#        else                                # COMMENT THIS
#            ifilt = 1:length(L1_f.lat);;    # COMMENT THIS
#        endif                               # COMMENT THIS

        L2.instr.(flds{i_instr}).data = L1_f.instr.(flds{i_instr}).data(ifilt,:);

        L2.instr.(flds{i_instr}).pitch = L1_f.instr.(flds{i_instr}).pitch(ifilt,:);
        L2.instr.(flds{i_instr}).roll = L1_f.instr.(flds{i_instr}).roll(ifilt,:);
        L2.instr.(flds{i_instr}).tilt = L1_f.instr.(flds{i_instr}).tilt(ifilt,:);
        L2.instr.(flds{i_instr}).(vza_nm) = L1_f.instr.(flds{i_instr}).(vza_nm)(ifilt,:);

        if isfield(L1_f.instr.(flds{i_instr}), 'int_time_sec');
            L2.instr.(flds{i_instr}).int_time_sec = L1_f.instr.(flds{i_instr}).int_time_sec(ifilt,:);
# blo: commented out as this variable is NOT in trios, so we don't want to add it if performing two filters
#        else 
#            L2.instr.(flds{i_instr}).int_time_sec = "";
        endif

        L2.instr.(flds{i_instr}).time = L1_f.instr.(flds{i_instr}).time(ifilt,:);

#        L2.instr.(flds{i_instr}) = hsas_cmp_stats(L2.instr.(flds{i_instr}));

    endfor



    flds = fieldnames(L2);

#    filter all data except for instrument data   
        for i_fld = 1:length(flds)

            if isstruct(L2.(flds{i_fld})) || size(L2.(flds{i_fld}),1)!=size(L1_f.gps.lat,1);
                continue
            endif

            L2.(flds{i_fld}) = L2.(flds{i_fld})(iLt,:);

        endfor

        
        L2.gps.time = L1_f.gps.time(ifilt,:);
        L2.gps.lat = L1_f.gps.lat(ifilt,:);
        L2.gps.lon = L1_f.gps.lon(ifilt,:);
        L2.gps.hdg = L1_f.gps.hdg(ifilt,:);





endfunction
