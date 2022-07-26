function hdr = hsas_rd_header(fn)

    for ifn = 1:length(fn)
#	disp(fn{ifn});
#	fflush(stdout);

        fid = fopen(fn{ifn});

            for irec = 1:10
                tmp = fgets(fid);
            endfor           
            
            tmp2 = strsplit(tmp, {" ","\t","\'"});  


            if length(tmp2)>2
          
                hdr.lat(ifn) = str2num(tmp2{2}) + str2num(tmp2{3})/60;

                if tmp2{4}=="S"
                    hdr.lat(ifn) = -hdr.lat(ifn);
                endif


                tmp = fgets(fid);
                tmp2 = strsplit(tmp, {" ","\t","\'"});
                
                hdr.lon(ifn) = str2num(tmp2{2}) + str2num(tmp2{3})/60;

                if tmp2{4}=="W"
                    hdr.lon(ifn) = -hdr.lon(ifn);
                endif

            else

                hdr.lat(ifn) = nan;
                hdr.lon(ifn) = nan;
            
            endif   


            for irec = 1:10
                tmp = fgets(fid);
            endfor           

            mesi = {"Jan",
                    "Feb",
                    "Mar",
                    "Apr",
                    "May",
                    "Jun",
                    "Jul",
                    "Aug",
                    "Sep",
                    "Oct",
                    "Nov",
                    "Dec"};

            tmp2 = strsplit(tmp, {" ","\t"}); 

   	    if length(tmp2)==6
 
              d = str2num(tmp2{4});
              m = find(strcmp(tmp2{3}, mesi));
              y = str2num(tmp2{6});
              tmp3 = strsplit(tmp2{5}, ":");
              H = str2num(tmp3{1});
              M = str2num(tmp3{2});
              S = str2num(tmp3{3});
    
              hdr.time(ifn) = datenum([y m d H M S]);
	
   	    else
		hdr.time(ifn) = nan;

	    endif


        fclose(fid);



    endfor








endfunction
