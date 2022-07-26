function tilt = cmp_tilt(roll, pitch)


    tilt = rad2deg(atan( sqrt(tan(deg2rad(roll)).^2 + tan(deg2rad(pitch)).^2) ) ) ;




endfunction