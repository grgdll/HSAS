function L1 = hsas_mk_L1_filtered(L0, max_tilt_accepted)

VERBOSE = false;



L1.wv = L0.wv;

iOK = find(L0.tilt < max_tilt_accepted & L0.instr.Li.data(:,1)>0 & L0.instr.Lt.data(:,1)>0 & L0.instr.Es.data(:,1)>0);

flds = fieldnames(L0);

for ifld = 2:length(flds)  # skip wv

    if isstruct(L0.(flds{ifld}))
        flds2 = fieldnames(L0.(flds{ifld}));

        for ifld2 = 1:length(flds2)

            if isstruct(L0.(flds{ifld}).(flds2{ifld2}))
                flds3 = fieldnames(L0.(flds{ifld}).(flds2{ifld2}));

                for ifld3 = 1:length(flds3)
                    
                    if size(L0.(flds{ifld}).(flds2{ifld2}).(flds3{ifld3}),1)==size(L0.(flds{ifld}).(flds2{ifld2}).data,1)
                        if VERBOSE, disp(["...L0." (flds{ifld}) "." (flds2{ifld2}) "." (flds3{ifld3})]); fflush(stdout); endif

                        L1.(flds{ifld}).(flds2{ifld2}).(flds3{ifld3}) = L0.(flds{ifld}).(flds2{ifld2}).(flds3{ifld3})(iOK,:);    
                    else
                        L1.(flds{ifld}).(flds2{ifld2}).(flds3{ifld3}) = L0.(flds{ifld}).(flds2{ifld2}).(flds3{ifld3});

                        if VERBOSE, disp(["..L0." (flds{ifld}) "." (flds2{ifld2}) "." (flds3{ifld3})]); fflush(stdout); endif
                    endif
    
                endfor

            else

                for ifld2 = 1:length(flds2)

                    if strcmp(flds{ifld}, "gps")
                        L1.(flds{ifld}).(flds2{ifld2}) = L0.(flds{ifld}).(flds2{ifld2})(iOK,:);  
                    elseif strcmp(flds{ifld}, "hdr")
                        L1.(flds{ifld}).(flds2{ifld2}) = L0.(flds{ifld}).(flds2{ifld2});
                    else
                        keyboard
                    endif
                    
                    if VERBOSE, disp([".L0." (flds{ifld}) "." (flds2{ifld2})]); fflush(stdout); endif
                endfor
    

            endif

        endfor

    else

        L1.(flds{ifld}) = L0.(flds{ifld})(iOK,:);

        if VERBOSE, disp(["L0." (flds{ifld})]); fflush(stdout); endif

    endif

endfor



endfunction
