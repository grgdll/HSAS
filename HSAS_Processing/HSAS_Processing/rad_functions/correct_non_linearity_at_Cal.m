function [out] = correct_non_linearity_at_Cal(rad_sn,sn_rad,sn,coeff,in,wv)
  disp("Non-linearity correction... Wait!")
  data = in;
  sn =str2num(sn);
  if sn == str2num(rad_sn.LI)
    coeff0=coeff.coeff_LI;
  endif
  if sn == str2num(rad_sn.LT)
    coeff0=coeff.coeff_LT;
  endif
  if sn == str2num(rad_sn.ES)
    coeff0=coeff.coeff_ES;
  endif

  alf = interp1(coeff0(:,1), coeff0(:,2), wv)';
  err = 1-alf.*data;
  data_corrected = data.*err;
  out = data_corrected;
endfunction