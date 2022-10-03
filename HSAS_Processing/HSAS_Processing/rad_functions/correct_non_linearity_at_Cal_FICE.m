function [out] = correct_non_linearity_at_Cal_FICE(rad_sn,sn_rad,sn,coeff,in,wv)
  # This is a modified version of correct_non_linearity_at_Cal for FICE 2022, accomodating new data format from tartu.
  # alpha is calculated within the function using S1 and S2 vectors and t1 and t2 integration times (refer to radcal_recipies.txt).

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


  # extract data from NL coeff matrix - refer to radcal_recipies.txt 
  wv_nl = coeff0(2:end,2); # nl subscript indicates wv vector is different length from nl matrix
  S1 = coeff0(2:end,7);
  S2 = coeff0(2:end,9);
  t1 = coeff0(1,7);
  t2 = coeff0(1,9);
  
  # compute alpha
  K = t1/(t2 + t1);
  S12 = (1 + K)*S1 - K*S2;
  alpha_nl = (S1 - S12)./S12.^2; 

  # previous steps duplicated from correct_non_linearity_at_Cal
  alpha = interp1(wv_nl, alpha_nl, wv,'extrap')'; # adjusts to length of cal file
  err = 1 - alpha.*data;
  data_corrected = data.*err;
  out = data_corrected;
  
  
endfunction
