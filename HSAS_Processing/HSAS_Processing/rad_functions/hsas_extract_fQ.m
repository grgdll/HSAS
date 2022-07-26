function L2 = hsas_extract_fQ(L2, alg)

    warning off
    addpath ../lib

    # load fQ table
        load("foQ.oct");
        L2.conf.foQ_table = "foQ.oct";
        
        if L2.conf.used_acs == true;
            L2 = hsas_read_ACS_Chl(L2); # L2.chl_acs added here
            L2.chl_estimated = L2.chl_acs;
        endif

        #if isfield(L2, "chl") & any(~isnan(L2.chl))# case for which we have chl data available
        if L2.conf.used_acs == true & any(~isnan(L2.chl_estimated)); # we use ACS Chl instead
            disp('ACS Chl used here')
            fflush(stdout);
            [L2.fQs0_fQs, L2.fQ0_fQs0, L2.fQs0, L2.fQs, L2.fQ0, L2.wv_ref] = hsas_fQ_ratio(L2.wv, L2.sza, L2.chl_estimated, L2.vza, L2.phi, foQ, PLOT=false, TEST=false);
        elseif isfield(L2, "chl_estimated");# case for which we have chl data available          
            disp("abbiamo chl_estimated");
            fflush(stdout);
            [L2.fQs0_fQs, L2.fQ0_fQs0, L2.fQs0, L2.fQs, L2.fQ0, L2.wv_ref] = hsas_fQ_ratio(L2.wv, L2.sza, L2.chl_estimated, L2.vza, L2.phi, foQ, PLOT=false, TEST=false);
        else
            # case for which we do not have chl data available and we iteratively estimate chl and fQ

               # initialize variables 
                    chl_old = 10;
                    i_iter = 1;

                #1) stima chla usando Lw/Ed
                    L2.chl_estimated = hsas_estimate_chl(L2.wv, L2.Rrs.data, alg);

                    innan = find(~isnan(L2.chl_estimated));  # indices of non-nan chl values

            # initialize results 
            L2.fQs0_fQs     = nan(size(L2.Rrs.data));
            L2.fQ0_fQs0     = nan(size(L2.Rrs.data));
            L2.fQs0         = nan(size(L2.Rrs.data));
            L2.fQs          = nan(size(L2.Rrs.data));
            L2.fQ0          = nan(size(L2.Rrs.data));
            L2.wv_ref       = nan(7,1);

            if any(innan)

                # set intial values of the i2p index (i2p=index to process)
                    i2p = find(abs(L2.chl_estimated-chl_old)>0.01 & ~isnan(L2.chl_estimated));
#                while all(abs(L2.chl_estimated(i2p)-chl_old)>0.01)
                while ~isempty(i2p)

  
                    #2) estrai fQis
                        [L2.fQs0_fQs(i2p,:), L2.fQ0_fQs0(i2p,:), L2.fQs0(i2p,:), L2.fQs(i2p,:), L2.fQ0(i2p,:), L2.wv_ref] = hsas_fQ_ratio(L2.wv, L2.sza(i2p,:), L2.chl_estimated(i2p,:), L2.vza(i2p,:), L2.phi(i2p,:), foQ, PLOT=false, TEST=false);

                    #3) compute LWN/Eo (exact normalized remote-sensing reflectance)
#                        LWN2E0.data = L2.Lw.data.*L2.R0_R.*L2.fQs0_fQs.*1./L2.instr.Es.avg.*L2.fQ0_fQs0;
                        LWN2E0.data = L2.Lw.data.*L2.R0_R.*L2.fQs0_fQs.*1./L2.instr.Es.data.*L2.fQ0_fQs0;

                    #4) updates values needed to check iteration
                        chl_old = L2.chl_estimated;
                        i_iter = i_iter + 1;
               
                    #5) stima chla usando LWN/Eo
                        L2.chl_estimated = hsas_estimate_chl(L2.wv, LWN2E0.data, alg);

                        i2p = find(abs(L2.chl_estimated-chl_old)>0.01 & ~isnan(L2.chl_estimated));


                            if i_iter>=9
                                disp("WARNING: max munber of interations reached")
                                fflush(stdout);
                                break
                            endif
                endwhile                
            else
                [L2.fQs0_fQs, L2.fQ0_fQs0, L2.fQs0, L2.fQs, L2.fQ0, L2.wv_ref, L2.chl_estimated] = deal([]);
            endif  
        endif
endfunction


