function out = rd_trios(fn)

disp(["reading file:  " fn]);
fflush(stdout);

    # get times
        fid = fopen(fn, 'r');
        tmp = fgetl(fid);
        fclose(fid);
        tmp_time = datenum(cell2mat(strtrim(strsplit(tmp, {"\#wv", " "}, true))(2:end-1)'), "HH:MM:SS");
        dbase    = strsplit(fn,"/"){end};
        dbase    = strsplit(dbase,"_"){end-2};
        out.time = datenum(dbase, "yyyy-mm-dd")+tmp_time-datenum(2017,1,1,0,0,0);
    # get data
        tmp = load(fn);
        out.wv = tmp(:,1);
        out.data= tmp(:,2:end)';

    if length(out.data(:,1))!=length(out.time)
        disp("\n\n---------Length of time array different from length of data: EXIT.\n\n")
        
        exit
    endif


disp("...done");
fflush(stdout);

endfunction
