function L2 = hsas_cmp_rho_2pars(L2)

VBS = false;
L2.rho_fitted = nan(length(L2.instr.Lt.data(:,1)),2);
#  L2.Lt0_fitted = nan(size(L2.instr.Lt.data(:,1)));

# define cost function
nir = find(L2.wv >=750 & L2.wv<=800); # [nm]
    
# fit rho
rho0   = 0.03;
Lt0ini = 0.5;
x0     = [rho0, Lt0ini];

record_length = length(L2.instr.Lt.data(:,1));

for id = 1:record_length;
   if VBS, disp(sprintf("%04u/%04u", id, length(L2.instr.Lt.data(:,1)))); fflush(stdout); endif
   L2.rho_fitted(id,:) = fminsearch(@(x) (sum(abs(  L2.instr.Lt.data(id,nir)-x(1)*L2.instr.Li.data(id,nir) -x(2) )) ), x0 );    
   rho0 = L2.rho_fitted(id,1);
   Lt0ini = L2.rho_fitted(id,2);
   x0 = [rho0, Lt0ini];
endfor

endfunction
