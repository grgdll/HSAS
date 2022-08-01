function L2 = hsas_extract_rho_AeronetOC(L2, ws_nm)


# set label to remember which table was used
    L2.conf.rho_table = "AeronetOC";
    L2.conf.rhoFOV = false;  # set this to false to remember when I used the whole FOV to compute rho


 ws = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0]';
vza = 40;
sza = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]';
phi = [90, 135]';


#case azimuth_offset of 90 
            rho(:,:,1)  =         [[0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256],# ;Wind speed of  0.0
                                        [0.0277,  0.0274,  0.0270,  0.0268,  0.0266,  0.0266,  0.0265,  0.0264],# ;Wind speed of  2.5
                                        [0.0357,  0.0309,  0.0285,  0.0281,  0.0278,  0.0277,  0.0276,  0.0274],# ;Wind speed of  5.0
                                        [0.0577,  0.0419,  0.0341,  0.0306,  0.0301,  0.0301,  0.0300,  0.0297],# ;Wind speed of  7.5
                                        [0.0802,  0.0597,  0.0440,  0.0344,  0.0333,  0.0328,  0.0326,  0.0323],# ;Wind speed of 10.0
                                        [0.1008,  0.0793,  0.0555,  0.0389,  0.0369,  0.0355,  0.0350,  0.0346],# ;Wind speed of 12.5
                                        [0.1180,  0.0963,  0.0678,  0.0437,  0.0406,  0.0380,  0.0371,  0.0363]]; #Wind speed of 15
#case azimuth_offset of 135 
            rho(:,:,2)  =         [[0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256,  0.0256],# ;Wind speed of  0.0
                                        [0.0270,  0.0268,  0.0267,  0.0267,  0.0268,  0.0268,  0.0266,  0.0264],# ;Wind speed of  2.5
                                        [0.0304,  0.0287,  0.0283,  0.0284,  0.0285,  0.0284,  0.0282,  0.0278],# ;Wind speed of  5.0
                                        [0.0406,  0.0321,  0.0304,  0.0304,  0.0307,  0.0306,  0.0303,  0.0296],# ;Wind speed of  7.5
                                        [0.0540,  0.0382,  0.0337,  0.0330,  0.0332,  0.0331,  0.0325,  0.0317],# ;Wind speed of 10.0
                                        [0.0679,  0.0471,  0.0378,  0.0361,  0.0356,  0.0355,  0.0347,  0.0337],# ;Wind speed of 12.5
                                        [0.0803,  0.0584,  0.0422,  0.0396,  0.0377,  0.0374,  0.0367,  0.0352]]; #Wind speed of 15


    ## find steps of rho axes
        s_ws = unique(diff(ws));
        s_sza = unique(diff(sza));


    # initialize output matrix
        L2.rho = nan(size(L2.gps.lat));

  


    # set ranges over which to interpolate
        i_ws = find_rng(ws, s_ws, L2.(ws_nm));
        i_sza = find_rng(sza, s_sza, L2.sza);



    for irec = 1:length(L2.gps.lat)    
        L2.rho(irec) = squeeze(interpn(ws(i_ws), sza(i_sza), squeeze(rho(i_ws,i_sza,1)), L2.(ws_nm)(irec), L2.sza(irec), 'linear'));
    endfor

L2.rho

    if any(isnan(L2.rho))
        disp("WARNING: some input variables out of the range used to generate rho tables");
        fflush(stdout);
    endif






endfunction


# find reduced ranges of values over which to search the LUT
function i_v = find_rng(v, s_v, v0)

   i_v = [];
   while length(i_v)<2   
       i_v = find(v>min(v0)-s_v & v<max(v0)+s_v);
       s_v = s_v*2;
   endwhile

       
   if min(v(i_v))>min(v0) | max(v(i_v))<max(v0)
       disp("PROBLEM")
       keyboard
   endif

endfunction

