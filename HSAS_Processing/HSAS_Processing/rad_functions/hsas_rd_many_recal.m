function out = hsas_rd_many_recal(sn, fn, VERBOSE)
# hsas = rd_hsas_many(sn, fn)
# 
# read many raw hsas files for one given instrument sn (ES, LI, LT)
 
    if VERBOSE, disp(["\nreading " sn{2} " files..."]); fflush(stdout); endif
    
    scaling_path='/data/datasets/cruise_data/active/AMT26/PML_optics/HEK_processing/code/hsas_recalibration/instrument_files/cal_files/cal_coeffs/';
    scaling_wavs=[1:1:1100];
    scale_vals  =ones([1,1100]);
    inter_wv    =[350:850];
    suffix='_coeffs_350850_1nm.txt';
    
    if strcmp(sn{2},'ES');
        scaling_file=[scaling_path 'ES' suffix];
    elseif strcmp(sn{2},'LI');
        scaling_file=[scaling_path 'LI' suffix];
    elseif strcmp(sn{2},'LT');
        scaling_file=[scaling_path,'LT' suffix];
    end
    fid = fopen(scaling_file,"r");
    file_vals = fscanf(fid,"%f", [1, Inf])';

    # pad out to 860
    scale_vals(inter_wv) = file_vals;
    scale_vals(isnan(scale_vals)) = 1.0;

    # remove nans
    scale_vals(isnan(scale_vals))=1.0;

    # find din
    idir = strfind(fn{1} , "/")(end);
    din = fn{1}(1:idir);

    out.data = [];
    out.time = [];
    out.int_time_sec = [];
    out.DARK_SAMPLE = [];
    out.DARK_AVE = [];
    out.TEMP_degC = [];
    out.FRAME_count = [];
    out.TIMER = [];
    out.CHECK_SUM = [];
    out.dark = [];
    out.recal_scaling = [];

    out.lat = [];
    out.lon = [];
    out.hdr_time = [];

    # load filenames of headers
    if strcmp(sn{2},"LT")
        fnhdr = glob( [din "????-???-??????.dat"] );
    endif
  
    for ifn = 1:length(fn)
        tmp = hsas_rd(fn{ifn});

        new_scale_vals = interp1(scaling_wavs,scale_vals,tmp.wv)';
        new_scale_vals = repmat(new_scale_vals,[size(tmp.(sn{2}))(1) 1]);
        out.data = [out.data; tmp.(sn{2})];
        out.recal_scaling = [out.recal_scaling; new_scale_vals];
        out.time = [out.time; tmp.time];
        out.int_time_sec = [out.int_time_sec; tmp.int_time_sec];
        out.DARK_SAMPLE = [out.DARK_SAMPLE,; tmp.DARK_SAMPLE];
        out.DARK_AVE = [out.DARK_AVE; tmp.DARK_AVE];
        out.TEMP_degC = [out.TEMP_degC; tmp.TEMP_degC];
        out.FRAME_count = [out.FRAME_count; tmp.FRAME_count];
        out.TIMER = [out.TIMER; tmp.TIMER];
        out.CHECK_SUM = [out.CHECK_SUM; tmp.CHECK_SUM];
        out.cal_file = tmp.cal_file;

        if fn{ifn}(:,end-8)=="D"        
            out.dark = [out.dark; ones(length(tmp.time),1)*true];
        else
            out.dark = [out.dark; zeros(length(tmp.time),1)*false];            
        endif

        # read header data (only during on sn = LT and not during dark counts)
        if strcmp(sn{2},"LT") 
            #find corresponding hdr file

            ifile = find(strcmp([fn{ifn}(1:end-12) ".dat"], fnhdr ));
            tmphdr = hsas_rd_header(fnhdr(ifile));
    
            out.hdr_time = [out.hdr_time; repmat(tmphdr.time, length(tmp.time),1)];
            out.lat = [out.lat; repmat(tmphdr.lat, length(tmp.time),1)];
            out.lon = [out.lon; repmat(tmphdr.lon, length(tmp.time),1)];

        endif

    endfor

        out.wv = tmp.wv;
        out.sv = tmp.sn;
endfunction
