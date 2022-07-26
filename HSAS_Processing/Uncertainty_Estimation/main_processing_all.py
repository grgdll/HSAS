import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from scipy import interpolate
import datetime
from pyephem_sunpath.sunpath import sunpos
from uncertainty_Rrs import uncertainty
import pandas as pd
from scipy.signal import savgol_filter
import pvlib
import pandas as pd
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import pickle


cruise = 'AMT26'
uncer = uncertainty(cruise)


#Folder including L1 and L2 data (Daily)
if cruise == 'AMT29':
    myrange = np.arange(287,326)
if cruise == 'AMT28':
    myrange = np.arange(268,302)
if cruise == 'AMT27':
    myrange = np.arange(266,308)
if cruise == 'AMT26':
    myrange = np.arange(267,307)

for i in myrange: #AMT26 267-307 #AMT27 266-308  #AMT28 268-302 #AMT29 287-326
    try:
        dir=str(i)

        #--Step1
        # Load L1, L2 data and some information
        uncer.load(dir)

        #-- Step2 detect clear sky
        uncer.sky_detect(is_figure=False)
        # --Step3
        # Uncertainty Calculation for L2 data
        uncer.perform_unc_summation()

        #---Check if save foler exist:
        save_filename = uncer.in_dir_L2.replace('L2','UncV4')+dir
        if not os.path.exists(save_filename):
            os.makedirs(save_filename)


        # ----Show results
        res = uncer.uncer_result
        
        lats = uncer.L2['L2'].gps.lat
        lons = uncer.L2['L2'].gps.lon
        dt = res['Hsas_dt']
        wl = res['wv']
        Rrs = res['Rrs']
        nLw = res['nLw']
        sig_Rrs = res['unc']
        sig_nLw = res['sig_nLw']
        air_temp = res['air_temp']
        sig_nLw = 100 * np.sqrt(np.array(sig_nLw)) / nLw
        day = i * 1.0
        clearsky = uncer.L2_clearsky
        viewphi = uncer.L2['L2'].phi#(uncer.L2['L2'].vaa)
        viewtheta = uncer.L2['L2'].vza
        suntheta = uncer.L2['L2'].sza
        windspeed = uncer.L2['L2'].true_wind_spd
        print(dir)
        print(save_filename)

        uncdata = {}
        uncdata['Rrs'] = (('datetime','wavelength'),Rrs)
        uncdata['sig_Rrs(%)'] = (('datetime','wavelength'),sig_Rrs)
        uncdata['nLw'] = (('datetime','wavelength'), nLw)
        uncdata['sig_nLw(%)'] = (('datetime','wavelength'),sig_nLw)
        datetime
        ds = xr.Dataset(
            uncdata,
            coords={
                'datetime': dt,
                'wavelength':wl,
                'air_temp':air_temp,
                'longitude': lons,
                'latitude': lats,
                'if_clearsky':clearsky,
            },
        )
        ds.to_netcdf(save_filename + '/AMT_Uncertainty_v4.0.nc', engine='netcdf4')
        # with open(r'uncdataAMT28.dat','wb') as output_file:
        #     pickle.dump(uncdata,output_file)
    except Exception as e:
        print(e)