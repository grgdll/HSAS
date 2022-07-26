function [v0 i_v] = hsas_find_rng(v, s_v, v0)
# find reduced ranges of values over which to search the LUT

    i_v = [];
    while length(i_v)<2   
        i_v = find(v>min(v0)-s_v & v<max(v0)+s_v);
        s_v = s_v*2;
    endwhile

        
#    if min(v(i_v))>min(v0) | max(v(i_v))<max(v0)
#        disp(sprintf("%s variable out of available range (%f): stopping", argn(1,:), mean(v0)))
#        exit
#    endif
        
    if min(v(i_v))>min(v0)
        disp(sprintf("%s variable out of available range (%f): setting it to min value", argn(1,:), min(v0)))
        v0(find(v0<min(v))) = min(v)*1.05;
    endif
        
    if max(v(i_v))<max(v0)
        disp(sprintf("%s variable out of available range (%f): setting it to max value", argn(1,:), max(v0)))
        v0(find(v0>max(v))) = max(v)*0.95;
    endif

endfunction


