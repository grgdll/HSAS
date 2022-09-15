function L2 = hsas_extract_rho_Mobley99(L2, ws_nm)
# function to compute rho using Mobley 1999 windspeed formula
# rho = 0.0256 + 0.00039 * windSpeed + 0.000034 * windSpeed^2


# set label to remember which method used
    L2.conf.rho_table = "Mobley_windspeed_formula";

#  from Mobley 1999
    L2.rho = 0.0256 + 0.00039*L2.appwind_spd + 0.000034*(L2.appwind_spd).^2


endfunction

