function L2 = hsas_R0_R(L2)
#
# function R0_R = hsas_R0_R(vza, wind_speed)
#
# compute ratio of gothic R at nadir to that at a given viewing zenith angle in air
# where theta is the viewing zenith angle in air {degrees], and wind_speed is in [m/s]
# 
warning off
addpath ../hsas.source/lib

    #### convert vza to theta_water
    # theta_water =  viewing zenith angle in water.   sin(vza) = n*sin(theta_water)   , where n=1.34
    # apply Snell Law    
        n = 1.34; # refractive index of seawater
        theta_water = rad2deg(asin(sin(deg2rad(L2.vza))/n));  # see above [degrees]
        

    ##### this table was downloaded in May 2016, but after converting d(:,2) to vza
    ##### I've plotted R vs vza and it seems that the values are still the wrong values that were corrected in Gordon 2005
    ##### Therefore, I am going to use R for a flat surface which is the closest to the true values
    
    d = load ("rgoth_Ebuchi.dat");
    L2.conf.gothR_table = "rgoth_Ebuchi.dat";    

    #  extract R(flat)

        iflat = find(d(:,1)<0);

        R_table = d(iflat,3);
        theta_water_table = d(iflat,2);            



    # interpolate R_table to theta values to be extracted

        L2.R = interp1(theta_water_table, R_table, theta_water);    

        L2.R0 = interp1(theta_water_table, R_table, theta_water*0);
    

    # compute 
        L2.R0_R = L2.R0./L2.R;




endfunction







