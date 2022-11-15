function chl = hsas_estimate_chl(wv, rrs, alg)

# tjor note: for FICE 2022, we use Adriatic Chl parameterization (specified by flag in proc_L2). This is based on 
# equation (2) in Berthon & Zibordi 2004.


    chl = nan(size(rrs,1),1);


    i490 = find(wv>=490,1);
    i555 = find(wv>=555,1);

#    if any( rrs(:,i490)<0 | rrs(:,i555)<0  )    
#        disp("negative rrs: can't estimate chl: exit");
#        chl = nan(size(rrs,1),1);
#        return
#    endif

    inan = find( rrs(:,i490)<0 | rrs(:,i555)<0  )  ;  
    rrs(inan,:) = nan;





  switch alg
    case { "K", "Black"}
      ialg = 1;
    case {"A", "Atlantic"}
      ialg = 2;
    case {"Adriatic"}
      ialg = 3;
    otherwise
      error ("invalid value");
  endswitch

#                   Black,      Atlantic,   Adriatic
chl_coeff =     [   -0.0067,      0.3190,     0.0910 
                    -2.6815,     -2.3360,    -2.6200 
                     1.2318,      0.8790,    -1.1480  #  tjor - this was previously +ve. It is -ve in B & Z 2004
                    -3.2713,     -0.1350,    -4.9490
                          0,     -0.0710,     0            ];
    

    RRS = log10(rrs(:,i490)./rrs(:,i555));


    chl = 10.^(chl_coeff(1) + chl_coeff(2,ialg)*RRS + chl_coeff(3,ialg)*RRS.^2 + chl_coeff(4,ialg)*RRS.^3 + chl_coeff(5,ialg)*RRS.^4);

    

endfunction




