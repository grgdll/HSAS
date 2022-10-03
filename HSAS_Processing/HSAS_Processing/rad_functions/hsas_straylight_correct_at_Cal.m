function data = hsas_straylight_correct_at_Cal(sensor_id,wl,sn_rad, D, data)
input_parameters_hsas;

# tjor-  matrix D is the straylight distribution matrix. For FICE, this has already been converted from the linear spread function in hsas calibrabate.m when the SL files are read

VBS = false;

I = eye(length(D),length(D));  
A = I + D;


disp(['** Straylight correction for ', sn_rad.(sensor_id), ' wait!'])

out = ones(size(data));
for i = 1:length(data(:,1))
 	correct_data = inv(A) * [data(i,:)]'; 
 	out(i,:) = correct_data(1:end);
end


data = out;


endfunction





