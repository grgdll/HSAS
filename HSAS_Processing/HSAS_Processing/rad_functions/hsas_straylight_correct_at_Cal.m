function data = hsas_straylight_correct_at_Cal(sensor_id,wl,sn_rad, D, data)
input_parameters_hsas;

VBS = false;

D = D(14:151,14:151); % hard coding step - 14 and 150 are the first and last wavelength bins of wl vector. bin 151 is for padding. This trunction step should be revisited when looking at class-based files.


I = eye(length(D),length(D)); % 

A = I + D;

disp(['** Straylight correction for ', sn_rad.(sensor_id), ' wait!'])

out = ones(size(data));
for i = 1:length(data(:,1))
 	correct_data = inv(A) * [data(i,:),0]'; # pads with a zero digit at the end
 	out(i,:) = correct_data(1:end-1);
end

data = out;

endfunction


