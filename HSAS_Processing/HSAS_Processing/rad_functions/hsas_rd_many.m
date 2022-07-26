function out = hsas_rd_many(sn, fn, VERBOSE)
# hsas = rd_hsas_many(sn, fn)
# 
# read many raw hsas files for one given instrument sn (ES, LI, LT)
 
    if VERBOSE, disp(["\nreading " sn{2} " files..."]); fflush(stdout); endif

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

    out.lat = [];
    out.lon = [];
    out.hdr_time = [];

    # load filenames of headers
    if strcmp(sn{2},"LT")
        fnhdr = glob( [din "????-???-??????.dat"] );
    endif
  


    for ifn = 1:length(fn)
		
		
		disp([fn{ifn}]);
		fflush(stdout);
		
        tmp = hsas_rd(fn{ifn});

        out.data = [out.data; tmp.(sn{2})];
        out.time = [out.time; tmp.time];
        out.int_time_sec = [out.int_time_sec; tmp.int_time_sec];
        out.DARK_SAMPLE = [out.DARK_SAMPLE,; tmp.DARK_SAMPLE];
        out.DARK_AVE = [out.DARK_AVE; tmp.DARK_AVE];
        out.TEMP_degC = [out.TEMP_degC; tmp.TEMP_degC];
        out.FRAME_count = [out.FRAME_count; tmp.FRAME_count];
        out.TIMER = [out.TIMER; tmp.TIMER];
        out.CHECK_SUM = [out.CHECK_SUM; tmp.CHECK_SUM];
        out.cal_file = tmp.cal_file;

		if strcmp(sn{2}, "ES") & strfind(strsplit(fn{ifn}, "/"){end}, "HED")       
		    out.dark = [out.dark; ones(length(tmp.time),1)*true];
		elseif (strcmp(sn{2}, "LT") | strcmp(sn{2}, "LI")) & strfind(strsplit(fn{ifn}, "/"){end}, "HLD")   
			out.dark = [out.dark; ones(length(tmp.time),1)*true];
		else
		    out.dark = [out.dark; zeros(length(tmp.time),1)*false];            
		endif

        # read header data (only during on sn = LT and not during dark counts)
        if strcmp(sn{2},"LT") 
            #find corresponding hdr file
            ifile = find(strcmp([fn{ifn}(1:end-12) ".dat"], fnhdr ));
				
		    if isempty(ifile)
				disp(['WARNING: missing header file ' [fn{ifn}(1:end-12) ".dat"]]);
				fflush(stdout);
				
	                tmphdr.time = nan;
	                tmphdr.lat = nan;
	                tmphdr.lon = nan;
		    else 
	            	tmphdr = hsas_rd_header(fnhdr(ifile));

			   	 if ~isstruct(tmphdr)
					keyboard
					disp(['WARNING: no data in header file ' fnhdr(ifile){1}]);
					fflush(stdout);
					tmphdr.time = nan;
					tmphdr.lat = nan;
					tmphdr.lon = nan;
				 endif
		    endif

            out.hdr_time = [out.hdr_time; repmat(tmphdr.time, length(tmp.time),1)];
            out.lat = [out.lat; repmat(tmphdr.lat, length(tmp.time),1)];
            out.lon = [out.lon; repmat(tmphdr.lon, length(tmp.time),1)];

        endif

    endfor

        out.wv = tmp.wv;
        out.sn = tmp.sn;




endfunction
