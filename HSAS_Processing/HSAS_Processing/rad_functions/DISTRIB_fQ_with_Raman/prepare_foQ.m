# prepare foQ matrix for octave
clear all

wv = [412 442 490 510 560 620 660]'; # wavelength [nm]
n_wv = length(wv);

sza = [0 15 30 45 60 75]'; # sun zenith angle
n_sza = length(sza);

vza = [1.078 3.411 6.289 9.278 12.300 15.330 18.370 21.410 24.450 27.500 30.540 33.590 36.640 39.690 42.730 45.780 48.830]';  # nadir_angles (in water?)
n_vza = length(vza);

vaa = [0 15 30 45 60 75 90 105 120 135 150 165 180]'; # azimuth_angles
n_vaa = length(vaa);

chl = [0.03 0.10 0.30 1 3 10]'; # chl concentration [mg/m3]
n_chl = length(chl);

# initialize foQ matrix

foQ = nan(n_wv, n_sza, n_vza, n_vaa, n_chl);


for i_wv = 1:n_wv
    for i_sza = 1:n_sza
        for i_chl = 1:n_chl
            fn =  glob(["/Users/gdal/Dropbox/Octave_functions/DISTRIB_fQ_with_Raman/D_foQ_pa/" sprintf("foQ_%03u\*_%02u_%05.2f", wv(i_wv), sza(i_sza), chl(i_chl)) ]){1} ;
            foQ(i_wv, i_sza, :, :, i_chl) = load(fn);
        endfor
    endfor
endfor


save -binary foQ.oct foQ


