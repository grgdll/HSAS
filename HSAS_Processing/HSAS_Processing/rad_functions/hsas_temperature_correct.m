function data = hsas_temperature_correct(data)

input_parameters_hsas;

VBS = false;

#-- number(#),wavelength(wl),correction factor(c), uncertianty(u)
v_es = dlmread([DIR_TEMPCORR,FN_TEMPCORR_ES],',',1,0);
v_lt = dlmread([DIR_TEMPCORR,FN_TEMPCORR_LT],',',1,0);
v_li = dlmread([DIR_TEMPCORR,FN_TEMPCORR_LI],',',1,0);

#---- interp to data(L0 or L1) wavelength
wv = data.wv;
c_es = interp1(v_es(:,2),v_es(:,3),wv);
c_lt = interp1(v_lt(:,2),v_lt(:,3),wv);
c_li = interp1(v_li(:,2),v_li(:,3),wv);
# u_es = interp1(v_es(:,2),v_es(:,4),wv);
# u_lt = interp1(v_lt(:,2),v_lt(:,4),wv);
# u_li = interp1(v_li(:,2),v_li(:,4),wv);

#-- read temperature 
##--*Note, for AMT cruises, temperature PCB in raw data are far incorrect as values are too large (e.g., 150 degree)*
##--*Thus,we used air temperautre instead*
#if isfield(data.surf_met,'air_temp1') & isfield(data.surf_met,'air_temp2')
#    temperature = (data.surf_met.air_temp1 + data.surf_met.air_temp2)/2;
#elseif isfield(data.surf_met,'air_temp1')
#    temperature = data.surf_met.air_temp1;
#elseif isfield(data.surf_met,'air_temp2')
#    temperature = data.surf_met.air_temp2;
#else
#    temperature = 25;
#endif


#---- read temperature
#---- We use pcb temperature here nstr.Lt.temp_pcb
if isfield(data.instr.Es,'temp_pcb')
  temp_Es = 0.474 * data.instr.Es.temp_pcb  - 48.539;
else
  temp_Es = 25;
endif
if isfield(data.instr.Lt,'temp_pcb')
  temp_Lt = 0.456 * data.instr.Lt.temp_pcb  - 45.438;
else
  temp_Lt = 25;
endif
if isfield(data.instr.Li,'temp_pcb')
  temp_Li = 0.462 * data.instr.Li.temp_pcb  - 45.934;
else
  temp_Li = 25;
endif

#----temperature correction
data.instr.Es.data = (1 - 0.001*c_es'.*(temp_Es - 21)).* data.instr.Es.data;
data.instr.Lt.data = (1 - 0.001*c_lt'.*(temp_Lt - 21)).* data.instr.Lt.data;
data.instr.Li.data = (1 - 0.001*c_li'.*(temp_Li - 21)).* data.instr.Li.data;
data.instr.Es.temp_pcb = temp_Es;
data.instr.Lt.temp_pcb = temp_Lt;
data.instr.Li.temp_pcb = temp_Li;
endfunction


