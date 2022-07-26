function [month,day] = jday2mmdd(year,doy)
	
    if is_leap_year(year)
        day_array = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        day_array = [31,28,31,30,31,30,31,31,30,31,30,31];
    endif
    sumday = 0;
    for i=1:12
        sumday = sumday + day_array(i);
        if sumday > doy
            month = i;
            day = doy - sumday + day_array(i);
            break
        endif
    endfor

endfunction
