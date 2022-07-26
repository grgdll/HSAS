function jday = jday(t)
	
    year = datevec(t)(:,1);
    month=datevec(t)(:,2);
    day = datevec(t)(:,3);
    
    if is_leap_year(year)
        day_array = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        day_array = [31,28,31,30,31,30,31,31,30,31,30,31];
    endif
    
    if i > 1
        jday = sum(day_array(1:month-1)) + day;
	else
        jday = day;
	endif
	
endfunction
