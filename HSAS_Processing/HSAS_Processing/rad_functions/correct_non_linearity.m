function [out] = correct_non_linearity(rad_sn,sn,coeff,in)
  disp("Non-linearity correction... Wait!")
  wv = in.raw_nodk.wv;
  data = in.raw_nodk.data;
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
  in.raw_nodk.data = data_corrected;
  out = in;
endfunction
