function [fQs0_fQs, fQ0_fQs0, fQs0, fQs, fQ0, wv] = hsas_fQ_ratio(wv0, sza0, chl0, vza0, phi0, foQ, PLOT, TEST);


VBS = false;


if ~exist("PLOT")
    PLOT = false;
endif

if ~exist("TEST")
    TEST = false;
endif



# initialize output matrix
     fQs = nan(length(vza0), length(wv0));
    fQs0 = nan(length(vza0), length(wv0));

     fQ = nan(length(vza0), length(wv0));
    fQ0 = nan(length(vza0), length(wv0));

    

#foQ = (n_wv, n_sza, n_vza, n_phi, n_chl);
# axes of foQ (in order)
     wv = [412.5 442.5 490 510 560 620 660]'; # wavelength [nm]
    sza = [0 15 30 45 60 75]'; # sun zenith angle
    vza = [1.078 3.411 6.289 9.278 12.300 15.330 18.370 21.410 24.450 27.500 30.540 33.590 36.640 39.690 42.730 45.780 48.830]';  # nadir_angles (in water?)
    phi = [0 15 30 45 60 75 90 105 120 135 150 165 180]'; # azimuth_angles
    chl = [0.03 0.10 0.30 1 3 10]'; # chl concentration [mg/m3]








## Marco interpola in quest'ordine
#phi, vza, sza, chl
    

if ~TEST
    #### convert theta_air to theta_water
    # theta_water =  viewing zenith angle in water.   sin(theta_air) = n*sin(theta_water)   , where n=1.34
    # apply Snell Law    
    n = 1.34; # refractive index of seawater
    vza_water0 = rad2deg(asin(sin(deg2rad(vza0))/n));  # see above [degrees]

else
    vza_water0 = vza0;  # see above [degrees]   USE THIS TO CHECK VALUES AGAINST Morel's table
    disp("WARNING: using input vza without refracting it into the water")

endif

# warn if wavelengths are outside of the tabulated range
#    if any(wv0>660 | wv0<412)
#        disp("Warning: wavelength range outside tabulated range => extrapolating")
#    endif
        
# warn if chl values are outside of the tabulated range
    if any(chl0>10 | chl0<0.03)
        disp("WARNING: chl range outside tabulated range => extrapolating")
    endif




# find steps of fQ axes
    s_sza = unique(diff(sza));
    s_phi = unique(diff(phi));

# set ranges over which to interpolate
    [sza0 i_sza] = hsas_find_rng(sza, s_sza, sza0);
    [phi0 i_phi] = hsas_find_rng(phi, s_phi, phi0);


iwv = find(wv0>=wv(1) & wv0<=wv(end));
# interpolate
    method = 'linear';

    for irec = 1:length(chl0)

        if VBS, disp( sprintf("%u/%u", irec,length(vza0)) ); fflush(stdout); endif

        if any(isnan([sza0(irec), vza_water0(irec), phi0(irec), log10(chl0)(irec)]))

            [sza0(irec), vza_water0(irec), phi0(irec), log10(chl0)(irec)]
            disp(["setting results to NaN at record " num2str(irec)]); fflush(stdout); 

            [fQs(irec, :), fQs0(irec, :), fQ0(irec, :)] = deal([nan]);
            continue
        endif


#[min(wv0(iwv)), min(sza(i_sza)),          min(vza), min(phi(i_phi)),   min(log10(chl));
#0,             sza0(irec),        vza_water0(irec),      phi0(irec), log10(chl0)(irec);
#max(wv0(iwv)), max(sza(i_sza)),           max(vza), max(phi(i_phi)), max(log10(chl))   ]

#keyboard
         fQs(irec, :) = squeeze(interpn(wv, sza(i_sza), vza, phi(i_phi), log10(chl), foQ(:,i_sza,:,i_phi,:), wv0, sza0(irec), vza_water0(irec), phi0(irec), log10(chl0)(irec), method, extrapval = nan));
        fQs0(irec, :) = squeeze(interpn(wv, sza(i_sza), vza, phi(i_phi), log10(chl), foQ(:,i_sza,:,i_phi,:), wv0, sza0(irec), 1.08,             phi0(irec), log10(chl0)(irec), method, extrapval = nan));
         fQ0(irec, :) = squeeze(interpn(wv, vza, phi(i_phi), log10(chl), squeeze(foQ(:,1,:,i_phi,:)),        wv0,             1.08,             phi0(irec), log10(chl0)(irec), method, extrapval = nan));

    endfor

# compute ratios
    fQs0_fQs = fQs0./fQs;

    fQ0_fQs0 = fQ0./fQs0;



if PLOT
    clf
    hold on
        plot(vza0, fQ)
#        plot(vza0, fQ*0+0.1, 'k')
#        plot(vza0, fQ*0+0.15, 'k')

#    legend(num2str(wv0'))
    ylim([0.05 0.2])
    grid on
endif




endfunction




