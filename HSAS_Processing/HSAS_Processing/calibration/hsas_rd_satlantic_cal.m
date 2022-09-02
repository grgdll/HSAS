function cal = hsas_rd_satlantic_cal(fn, used_pixels)
keyboard
    if (nargin == 1)
        used_pixels = 1:255;
    endif


    fid = fopen(fn, 'r');
    
    tmp = '         ';
    
# skip header
    while strcmp(tmp(1:6),"# Date")==0
        tmp = fgets(fid);
    endwhile
    fgets(fid);

#read date of cal, instrument type, integration time
    tmp = fgets(fid);   
# # 2017-05-05-18-53-51 |i_ansko  |1.1        |A   |ES   |256.0    |
    tmp = strsplit(tmp,{' ', '\t', '|'});
    
    cal.date = datenum(tmp{2}, "yyyy-mm-dd-HH-MM-SS");
    cal.instru = tmp{6};
    cal.int_time = num2str(tmp{5}); # [milliseconds]
    
 # skip more header
    tmp = ' ';
    while isempty(strfind(tmp,"INSTRUMENT"))
        tmp = fgets(fid);
    endwhile
    tmp = strsplit(fgets(fid)); 
    
    cal.sn = tmp{2};
    
    
    
    
# skip more header
    tmp = ' ';
    while isempty(strfind(tmp,"# Spectrum Data"))
        tmp = fgets(fid);
    endwhile
    fgets(fid);
    
    
  # Number of Dark Samples  
    
 # read calibration file
 
# ES 305.66 'uW/cm^2/nm' 2 BU 1 OPTIC3
# 997.091	0.00000000E+000	1.000	0.256
# 
    clear tmp
    tmp = fgets(fid);
    iwv = 1;
    while isempty(strfind(tmp,"# Number of Dark Samples"))
        tmp = strsplit(tmp); # read wavelength
        tmp2 = strsplit(fgets(fid));  # read calibration
        fgets(fid); # skip empty line
		
        cal.wv(iwv) = str2num(tmp{2});
        cal.offset(iwv) = str2num(tmp2{1});
        cal.gain(iwv) = str2num(tmp2{2});
        cal.int_time_wv(iwv) = str2num(tmp2{4}); # [seconds]
        
        # convert tmp back into one string
#         tmp = strrep(cell2mat(strcat(tmp, '-')), '-', ' ')
        tmp = fgets(fid);
        
#         if strfind(tmp, '1142.23 ')
#         keyboard
#         endif
        
        iwv = iwv+1;
        
        
    endwhile
   

    
    
    fclose(fid);


    cal.wv = cal.wv(used_pixels);
    cal.offset = cal.offset(used_pixels);
    cal.gain = cal.gain(used_pixels);
    cal.int_time_wv = cal.int_time_wv(used_pixels);






endfunction



