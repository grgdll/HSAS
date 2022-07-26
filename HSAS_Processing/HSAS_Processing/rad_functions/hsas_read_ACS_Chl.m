function L2 = hsas_read_ACS_Chl(L2)

input_parameters_hsas;

data = importdata(FN_ACS);
chl = [];
datetime_num = [];
for i=2:length(data); 
  strs = strsplit(data{i},',');
  chl = [chl,str2num(strs{20})];
  #2016-09-23T07:19:07Z
  strarray = sscanf(strs{3},'%d-%d-%dT%d:%d:%dZ');
  datetime_num = [datetime_num, datenum(strarray(1),strarray(2),strarray(3),strarray(4),strarray(5),strarray(6))];
endfor

L2.chl_acs = interp1(datetime_num,chl,L2.time);

endfunction