function L2 = hsas_extract_F0(L2)

#    F0 = load("./lib/Thuillier_F0.txt");
##    F0 = load("./lib/E0_giorgio.txt");

#    wv = F0(:,1);
#    F0smooth = slidefun(@mean, 11, F0(:,2));  # Marco usa una moving average che da' lo stesso risultato di questa funzione

#    F0 = interp1(wv, F0smooth, L2.wv)';
#    save -binary ./lib/E0_giorgio.oct F0

#    keyboard

    
    warning off
    addpath ../lib

    load("./rad_functions/E0_giorgio.oct");
    L2.conf.E0_table = "E0_giorgio.oct";
    
    #L2.F0 = F0*10;  # [mW/m^2/nm]  CONVERT TO TRIOS UNITS 
    L2.F0 = F0; # [uw/cm^2/nm] Changed by Junfang
endfunction
