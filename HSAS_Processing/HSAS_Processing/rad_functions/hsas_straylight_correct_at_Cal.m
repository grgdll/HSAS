function data = hsas_straylight_correct_at_Cal(sensor_id,wl,sn_rad, D, data)
input_parameters_hsas;

VBS = false;

I = eye(length(wv),length(wv));

#----Correction for Es
A = I +D;
disp(['** Straylight correction for ', sn_rad.(sensor_id), ' wait!'])


out = ones(size(data));
for i = 1:length(data(:,1))
 correct_data = inv(A) * [data(i,:),0]';
 out(i,:) = correct_data(1:end-1);
end

data = out;

endfunction


