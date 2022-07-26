#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:01:49 2022

@author: tjor

Python implementation of HSAS radiometic post-processing & data formatting for 
FICE 2022, Adriatic

Script includes:

(i) Loading of L1 (Es, Li, Lt) and L2 (Rrs) output .mat files and end-end 
uncertainty estimates

(ii) Spectral downsampling of Li, Lt, Es, Rrs using OLCI spectral response function

(iii) Writing output spectra (downsampled to OLCI bands) to csv file (for each station)
- Includes Rrs, Es, Li, Lt, & metadata for each timestamp 
- Metadata includes: rho, wind-speed, `cloudiness index: pi*Li(400)/Es(400), 
mask for glint (Lt-based) QC, mask for tilt QC, mask for,

(iv) Station-averaging (mean, median, STD) of radiometric and wind data

(v) Writing station summary to .csv file based on FICE_2022_AAOT template
- Includes station avergaes for Rrs, Es, Lt, Lt, Wind & core metadata required

"""
#"""


import os
import sys
import csv
import netCDF4 as nc
import pandas as pd

import numpy as np
import glob   
import pandas as pd
import ephem

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import math
import matplotlib.dates as mdates
import datetime as dt
from scipy.interpolate import interp1d
import scipy.io as io
import datetime

 #  functions to unpack srf and octave files in np array format
def unpack_srf(dir_srf):
    '''Function to unpack OLCI SRF in np array format'''
   
    srf_data = nc.Dataset(dir_srf + '/S3A_OL_SRF_20160713_mean_rsr.nc4') # load OLCI SRF
    srf_bands =  srf_data['nominal_centre_wavelength'][0:19]  # band centres up to 900 nm
    # band_width = srf_data['bandwidth_fwhm'][0:16] - FWHM of SRF - not needed
    srf_wv = np.array(srf_data['mean_spectral_response_function_wavelength'])
    srf = np.array(srf_data['mean_spectral_response_function'])

    return srf, srf_wv, srf_bands


def unpack_L1(fn_L1):
    ''' unpacks relevant L1 data in np array format: L1['L1'].__dict__.keys()
    can be used to show fields - needs to be applied at each struct level'''
    
    L1 = io.loadmat(fn_L1[0], squeeze_me = True, struct_as_record = False) 
    
    time = np.array(numeric_to_UTC_time(L1['L1'].gps.time)) 
    windspeed = L1['L1'].appwind_spd
    phi = L1['L1'].phi
 
    Es = L1['L1'].instr.Es.data # hyperspectral data
    Lt = L1['L1'].instr.Lt.data
    Li = L1['L1'].instr.Li.data
    
    if 'int_time_sec' in L1['L1'].instr.Es.__dict__.keys():
        Es_int_t =  L1['L1'].instr.Es.int_time_sec[0] # integration time seconds - assumes constant for each station
        Lt_int_t =  L1['L1'].instr.Lt.int_time_sec[0]
        Li_int_t = L1['L1'].instr.Li.int_time_sec[0]
    else:
        Es_int_t =  np.nan
        Lt_int_t =  np.nan
        Li_int_t = np.nan
    
    wv = L1['L1'].wv 

    return time, windspeed, phi, Es, Lt, Li, Es_int_t, Lt_int_t, Li_int_t, wv


def unpack_L2(fn_L2, time_L1, wv):
    ''' unpacks relevant L2 data in np array format. Adjustment is made so that 
    size of output matches L1 data. Timestamps that do not pass quality control are NaN padded'''
  
    L2 = io.loadmat(fn_L2[0], squeeze_me = True, struct_as_record = False) 
    
    # convert L2 data to length of L1 data, with nan-padding added where no data
    qc_mask = np.zeros(len(time_L1)) # 1 == data, 0 === no data
    Rrs = np.nan*np.ones([len(time_L1), len(wv)])
    exLwn = np.nan*np.ones([len(time_L1),len(wv)])
    rho = np.nan*np.ones(len(time_L1))
    phi = np.nan*np.ones(len(time_L1))
    
    # load files in L2 format- L2 subscript indicates matrices are undersize
    if 'Rrs' in L2['L2'].__dict__.keys():
        Rrs_L2 = L2['L2'].Rrs.data 
        rho_L2 = L2['L2'].rho 
        exLwn_L2 = L2['L2'].exLwn.data
        time_L2 = np.array(numeric_to_UTC_time(L2['L2'].gps.time)) 
        phi_L2 = L2['L2'].phi
        
        # indicies in L1 time that match L2 time (i.e. pass overall QC)
        L1_matches = np.intersect1d(time_L1, time_L2, return_indices = True)[1]
        for i in range(len(L1_matches)): # fill L1-size stuctures 
            qc_mask[int(L1_matches[i])] = 1 
            Rrs[int(L1_matches[i]),:] = Rrs_L2[i,:] 
            exLwn[int(L1_matches[i]),:] = exLwn_L2[i,:] 
            rho[int(L1_matches[i])] = rho_L2[i]
            phi[int(L1_matches[i])] = phi_L2[i]
            
    return qc_mask, Rrs, exLwn, rho, phi


def unpack_L2_nomask(fn_L2):
    ''' Unpacks relevant L2 data in np array format. No mask - just includes 
    timestamps where Rrs retreivals were successful'''
  
    L2 = io.loadmat(fn_L2[0], squeeze_me = True, struct_as_record = False) 
    
    time = np.array(numeric_to_UTC_time(L2['L2'].gps.time)) 
    wv = L2['L2'].wv 
    
    # (Ir)radiance
    Es = 10*L2['L2'].instr.Es.data   # Factor 10 coverts from mw cm^-2 um^-1 to mW m^-2 nm^-1 
    Li = 10*L2['L2'].instr.Li.data 
    Lt = 10*L2['L2'].instr.Lt.data 
    
    if 'int_time_sec' in L2['L2'].instr.Es.__dict__.keys():
        Es_int_t = L2['L2'].instr.Es.int_time_sec[0] # integration time seconds - assumes constant for each station (which it is)
        Lt_int_t = L2['L2'].instr.Lt.int_time_sec[0]
        Li_int_t = L2['L2'].instr.Li.int_time_sec[0]
    else:
        Es_int_t = np.nan
        Lt_int_t = np.nan
        Li_int_t = np.nan
     
    # Rrs
    Rrs = L2['L2'].Rrs.data 
    exLwn = 10*L2['L2'].exLwn.data # Factor 10 coverts from mw cm^-2 um^-1 to mW m^-2 nm^-1 

    windspeed = L2['L2'].appwind_spd
    rho = L2['L2'].rho   
        
    return time, wv, Es, Li, Lt, Es_int_t, Lt_int_t, Li_int_t, Rrs, exLwn, windspeed, rho


def numeric_to_UTC_time(time_num):
     '''Function to convert numeric time (octave output for Hsas spectra) to UTC 
     - based on Junfangs' uncertainty code'''
   
     time_UTC = []
     for i in range(len(time_num)):
         time_UTC_i = (datetime.datetime.fromordinal(int(time_num[i])) + 
                        datetime.timedelta(seconds = 86400*(time_num[i] - float(int(time_num[i])))) - datetime.timedelta(days=366))
         
         time_UTC.append(time_UTC_i)
     
     return time_UTC


def hyperspec_to_OLCIbands(S, time, wv, srf, srf_bands, srf_wv):
    ''' Function to downsample hyperspectral data to OLCI bands: input spectra 
   in dimenions of time * wavelength'''
   
    D = np.nan*np.zeros([len(time),len(srf_bands)]) # matrix format for down-sampled spectra - time * spectral band

    # Loop for spectral down-sampling (outputs np array D)
    for j in range(len(S)): # loop over timestamps
        if np.sum(np.isnan(S[j,:])) < len(S[j,:]): # checks spectra exists ()
            wv_j = wv[np.isnan(S[j,:]) == 0] # remove wl bins with nan padding (rq. for interpolation functions)
            S_j = S[j, np.isnan(S[j,:]) == 0]    
            for k in range(len(srf_bands)): # loop over spectral bands
                first_wv = wv[int(np.where(np.isnan(S[j,:]) == 0)[0][0])]  # first and last wavelengths
                last_wv = wv[int(np.where(np.isnan(S[j,:]) == 0)[0][-1])]
                if np.min(srf_wv[k]) > first_wv and np.max(srf_wv[k]) < last_wv: # tests data exists in band
                    interp_funct = interp1d(wv_j, S_j, kind = 'cubic') # interpolate to same wv interval as OLCI
                    S_j_interp = interp_funct(srf_wv[k]) # oversample on same wv range as OLCI SRF 
                    S_j_k = np.sum(S_j_interp*srf[k])/np.sum(srf[k]) # convolution of kth band for jth timestamp with OLCI SRF
                    D[j,k] =  S_j_k # fill data matrix
    
    # Loop for pandas dataframe format
    D_df = pd.DataFrame(index = time) # pandas data frame format for down-sampled spectra
    for k in range(len(srf_bands)):
           D_df[str(srf_bands[k])] = D[:,k] 
                
    return D, D_df


def multi_hyper_plot(S, wv, D, srf_bands, ylab):
    ''' test plot for downsampling'''
   
    plt.figure(figsize=(10,10))
    plt.scatter(srf_bands,D[1].T[1], label = 'OLCI band centres')  
    plt.plot(wv,S[1,:], label = 'Hyperspectral', color = 'red')
    plt.legend()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel(ylab)            
    
    return 


def station_averages(D):
    ''' mean, median ,std'''
   
    D_mean = np.nanmean(D,axis=0)
    D_med = np.nanmedian(D,axis=0)
    D_std = np.nanstd(D,axis=0)        
    
    return D_mean, D_med, D_std


def station_averages_dataframe(D, spec_id, station_i, time_L1, windspeed, phi, srf_bands):
    '''computes averages in data frame format including station meta data '''
          
    summary_df = pd.DataFrame(index=None) 
    summary_df['spectrum'] = [spec_id + '_median', spec_id + '_standard_deviation', spec_id + '_mean', spec_id + '_uncertainty'] 
    summary_df['station'] = [station_i, station_i, station_i, station_i]
    summary_df['water body'] = ['Adriatic Sea', 'Adriatic Sea', 'Adriatic Sea', 'Adriatic Sea']
    summary_df['station start time'] = [str(time_L1[0])[0:19], str(time_L1[0])[0:19], str(time_L1[0])[0:19], str(time_L1[0])[0:19]]
    summary_df['wind speed [m/s] (median, std, mean)'] = [str(np.nanmedian(windspeed)), str(np.nanstd(windspeed)), str(np.nanmean(windspeed)), 'NaN']
    summary_df['phi'] = [phi[0], phi[0], phi[0], phi[0]]  
    
    # averages and std and uncertainty  of downsampled spectra in each band (np.array format)
    D_summary = np.zeros([4,len(srf_bands)])
    D_summary[0,:] = np.nanmedian(D,axis=0)
    D_summary[1,:] = np.nanstd(D,axis=0)
    D_summary[2,:] = np.nanmean(D,axis=0)
    D_summary[3,:] = np.nan*np.ones(len(srf_bands))  # temporary fix for uncertainty

    # fill data frame with np.array info
    for k in range(len(srf_bands)):
        summary_df[str(srf_bands[k])] = D_summary[:,k] 
  
    return summary_df


def write_spectra(D_df,sensor_info, fpath):
    ''' writes downsampled spectra to csvfile'''
    
    # create headers using sensor info dictionary
    with open(fpath, 'w', encoding = 'UTF8') as f:
        writer = csv.writer(f)
        for i in range(len(sensor_info)):
            row_i = str('# '  + list(sensor_info.keys())[i] + ': ' + str(list(sensor_info.values())[i]))
            print(row_i)
            writer.writerow([''.join(row_i)])
    
    # append spectra dataframe to csvfile: dim = time * OLCI bands
    D_df.to_csv(fpath, na_rep ='NaN', mode = 'a')
    
    return 

       
def write_station_summary(summary_path, time, Es_av, Lt_av, Li_av, Rrs_av, exLwn_av):
    ''' Function to write station summary in csv format'''
    
    # fill up headers
    with open(summary_path, 'w', encoding='UTF8') as f:
        writer = csv.writer(f)           
        row = str('#file_type: spectrums')
        writer.writerow([''.join(row)])
        row = str('#file_version: 1.8.0')
        writer.writerow([''.join(row)])
        row = str('#file_created: ' + str(dt.datetime.now())[0:19])
        writer.writerow([''.join(row)])
        row = str('#data_type: HyperSAS_above')
        writer.writerow([''.join(row)])
        row = str('#is_grouped_by_station: true')
        writer.writerow([''.join(row)])
        row = str('#correction: none')
        writer.writerow([''.join(row)])
        row = str('#wavelengths: S3A_OLCI')
        writer.writerow([''.join(row)])
        row = str('#time_period: [{"start:"' + str(datetime.date(2022, 7, 13)) + ',"end:"' + str(datetime.date(2022, 7, 21)) + '}]')
        # row = str('#time_period: [{"start:"' + str(time[0])[0:19] + ',"end:"' + str(time[-1])[0:19] + '}]')
        writer.writerow([''.join(row)])
        row = str('#reflectance_calculation equation: M99')
        writer.writerow([''.join(row)])
        row = str('units: {"reflectance":["Wavelength [nm]","Reflectance [sr^-1],"es":["Wavelength [nm]","Irradiance [mW m^-2 nm^-1]"],"lt":["Wavelength [nm]","Radiance [mW m^-2 nm^-1 sr^-1]"],"li":["Wavelength [nm]","Radiance [mW m^2 nm^-1 sr^-1)]", "exLnw":["Wavelength [nm]","Radiance [mW m^2 nm^-1 sr^-1)]"]}')
        writer.writerow([''.join(row)])
    
    # append station average dataframes for each sectrum to csv file -
    for i in range(len(Rrs_av)):
        Rrs_av[i].to_csv(summary_path, na_rep ='NaN', index = False, mode = 'a')
        Es_av[i].to_csv(summary_path, na_rep ='NaN', index = False, mode = 'a')
        Lt_av[i].to_csv(summary_path, na_rep ='NaN', index = False, mode = 'a')
        Li_av[i].to_csv(summary_path, na_rep ='NaN', index = False, mode = 'a')
        exLwn_av[i].to_csv(summary_path, na_rep ='NaN', index = False, mode = 'a')
        
    return


def append_to_summary(dict_list):
     ''' Function to append to overall summary'''
    
     # cloudiness index - pi*Li(400)/Es(400)
     env_index = np.pi*np.array(Li_OLCI[1]['400.0'])/np.array(Es_OLCI[1]['400.0']) 
     env_index_mean = np.nanmean(env_index)
     # env_index_std =  np.nanstd(env_index)
     
     N_rrs = np.sum(~np.isnan(np.array(Rrs_OLCI[1]['400.0']))) # no. of rrs retrievals
     # P_rrs = 100*N_rrs/len(Rrs_OLCI[1]['400.0'])     # % of rrs retreivals
     
     CV_rrs = np.nanstd(Rrs_OLCI[1],axis=0)/np.nanmean(Rrs_OLCI[1],axis=0)
     CV_rrs_band_av = 100*np.nanmean(CV_rrs[0:6])  # % of rrs retreivals
     
     windspeed_mean = np.nanmean(windspeed)
     # windspeed_std = np.nanstd(windspeed)
     
     # select data based on environmental conditions
     if  CV_rrs_band_av < 2 and env_index_mean < 0.4: # good data
         env_mask = 'Valid'
     else:                                            # poor data
         env_mask = 'Not Valid'
    
     new_row = {'station':  str(station_index[i]) + '_' + str(stations[i]), 'start time [UTC]': str(time[0]),
                'end time [UTC]': str(time[-1]), 'mask':  env_mask, 'env_index_mean': env_index_mean, 
                 'number_Rrs': N_rrs, 'CV_rrs_band_av': CV_rrs_band_av,
                 'windspeed_mean': windspeed_mean}
 
     dict_list.append(new_row)
    
     return dict_list
 
    
def conditions_summaryplots(df_overall):

     fig, ax = plt.subplots()
     plt.figure(figsize=(30, 20))
     plt.rc('font', size=20)   
     plt.title('Dependence of Rrs variability on cloudiness (station average values)')
     colors = cm.jet(np.linspace(0, 1, len(df_overall)))
     for i in range(len(df_overall)):
         plt.scatter(df_overall['CV_rrs_band_av'][i], df_overall['env_index_mean'][i],c=colors[i], cmap='jet',label = df_overall['station'][i])
     plt.legend(labels = df_overall['station'], ncol=2,fontsize=18)
     plt.xlabel('CV[Rrs]: mean for OLCI bands 1-6 [%]')
     plt.ylabel('$\pi$Li(400)/Es(400): (cloudiness index)')
     plt.ylim(0,1)
     plt.xlim(0,10)
     plt.savefig(dir_write + '/FICE_cloudiness_CVRrs.png')
     
     plt.figure(figsize = (30, 20))
     plt.rc('font', size=20)   
     plt.title('Dependence of Rrs variability on windspeed (station average values)')
     colors = cm.jet(np.linspace(0, 1, len(df_overall)))
     for i in range(len(df_overall)):
         plt.scatter(df_overall['CV_rrs_band_av'][i], df_overall['windspeed_mean'][i],c=colors[i], cmap='jet',label = df_overall['station'][i])
     plt.legend(labels = df_overall['station'], ncol=2, fontsize=18)
     plt.xlabel('CV[Rrs]: mean for OLCI bands 1-6 [%]')
     plt.ylabel('Windspeed [m/s]')
     plt.ylim(0,6)
     plt.xlim(0,10)   
     plt.savefig(dir_write + '/FICE_windspeed_CVRrs.png')
 
     return
    
    
    
if __name__ == '__main__':
    
    # subdirectories to read data
    dir_main = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/'
    # dir_procrun = 'Nocorrections_FICEcalfiles/FICE_secondsubmission_quadratic_windspeeds/' # 'processing run'
    dir_procrun = 'Nocorrections_FICEcalfiles/FICE_sixthsubmission_NIR_V2/'
    
    dir_cal = dir_main + 'Processed_TPcorrection/' + dir_procrun + 'Calibrated/'
    dir_L0 = dir_main + 'Processed_TPcorrection/' +  dir_procrun + 'L0/'
    dir_L1 = dir_main + 'Processed_TPcorrection/' +  dir_procrun + 'L1/'
    dir_L2 = dir_main + 'Processed_TPcorrection/' +  dir_procrun + 'L2/'
    dir_unc = dir_main + 'TO_ADD' 
    dir_srf = dir_main + '/Source/HSAS/HSAS_Processing/FICE_Formatting' # cuntains OLCI SRF and template
    
    # subdirectories to write data
    dir_write = dir_main +'Processed_TPcorrection/' + dir_procrun  + 'FICE_format/'
    if os.path.exists(dir_write) == False:
       os.mkdir(dir_write)
    
    # initialize station sub directories
    stations = sorted(os.listdir(dir_L0))
    station_index = np.arange(0,78,1)
    station_index[70:-1] = np.arange(71,78,1) # raw data from station 70 _134000 is incomplete - cant process
    for i in range(len(stations)):
      os.mkdir(dir_write + str(station_index[i]) + '_' + str(stations[i]))
      
    # load OLCI srf
    srf, srf_wv, srf_bands = unpack_srf(dir_srf)

    Es_av_summary = []
    Lt_av_summary = []
    Li_av_summary = []
    Rrs_av_summary = []
    exLwn_av_summary = []

    dict_list = []   # list to append summaries
    # for i in range(1): # process each station in sequence
    for i in range(len(stations)): # process each station in sequence
        # access  filenames (hsas data structures) of ith station - cal and L0 currently not used
        # fn_cal = glob.glob(dir_cal + stations[i]  + '/*dat*')
        fn_L0 = glob.glob(dir_L0 + stations[i]  + '/*mat*')        
        fn_L1 = glob.glob(dir_L1 + stations[i]  + '/*mat*') #
        fn_L2 = glob.glob(dir_L2 + stations[i]  + '/*mat*')
        
        # fig_L2 = glob.glob(dir_L2 + stations[i]  + '/plot/*png') # optional scrip - re-name plots
        # os.rename(fig_L2[0], dir_L2 + stations[i]  + '/' + str(station_index[i]) + '_HSAS_Level2summaryplots.png')
        
        # unpack variables from L1 and L2 data data stuctures in np array format
        time_L1, windspeed_L1, phi_L1, Es_L1, Lt_L1, Li_L1, Es_int_time_L1, Lt_int_time_L1, Li_int_time_L1, wv = unpack_L1(fn_L1)
        # qc_mask, Rrs, exLwn, rho, phi = unpack_L2(fn_L2, time, wv) # L2 data is just where rrs passes QC - it is nan-padded to match length of L1 data
        time, wv, Es, Li, Lt, Es_int_t, Lt_int_t, Li_int_t, Rrs, exLwn, windspeed, rho =  unpack_L2_nomask(fn_L2)

        # spectral downsampling to OLCI: putput element 0 is np array, 1 is dataframe format
        Es_OLCI = hyperspec_to_OLCIbands(Es, time, wv, srf, srf_bands, srf_wv)  
        Lt_OLCI = hyperspec_to_OLCIbands(Lt, time, wv, srf, srf_bands, srf_wv)
        Li_OLCI = hyperspec_to_OLCIbands(Li, time, wv, srf, srf_bands, srf_wv)
        Rrs_OLCI = hyperspec_to_OLCIbands(Rrs, time, wv, srf, srf_bands, srf_wv) 
        exLwn_OLCI = hyperspec_to_OLCIbands(exLwn, time, wv, srf, srf_bands, srf_wv) 
        # multi_hyper_plot(Es, wv, Es_OLCI, srf_bands, 'Es [mW cm$^{2}$ $\mu$m^{1}$')
        

        # collate sensor info and write downsampled spectra to csv files
        wv_string = '(' + str(wv[0]) + ', ' + str(wv[-1]) + ', ' + str(wv[1]-wv[0]) +')'
        Es_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '2027A', 
                   'Integration time': str(Es_int_t) + ' [sec]', 
                   'Wavelength scale (max, min, spacing: all post interpolation)':  wv_string,
                   'Units':  ' [mW m^-2 nm^-1]'}
                
        Lt_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '464', 
                   'Integration time': str(Lt_int_t) + ' [sec]', 
                   'Wavelength scale (max, min, spacing: all post interpolation)': wv_string,
                   'Units':  ' [mW m^-2 nm^-1 sr^-1]'}
                 
        Li_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '2054A', 
                   'Integration time': str(Li_int_t) + ' [sec]', 
                   'Wavelength scale (max, min, spacing: all post interpolation)': wv_string,
                   'Units':  ' [mW m^-2 nm^-1 sr^-1]'}
                       
        Rrs_info = {'Make': 'Seabird', 'Model': 'HypserSAS',
                    'Wavelength scale (max, min, spacing: all post interpolation)': wv_string,
                    'Units':  ' [sr^-1]'}
                        
        exLwn_info = {'Make': 'Seabird', 'Model': 'HypserSAS',
                      'Wavelength scale (max, min, spacing: all post interpolation)': wv_string,
                      'Units':  ' [mW m^-2 nm^-1 sr^-1]'}
                        
        Es_path = dir_write + str(station_index[i]) + '_' + str(stations[i]) + '/' + 'Es_' + str(station_index[i]) + '.csv'
        Lt_path = dir_write + str(station_index[i]) + '_' + str(stations[i]) + '/' + 'Lt_' + str(station_index[i]) + '.csv'
        Li_path = dir_write + str(station_index[i]) + '_' + str(stations[i]) + '/' + 'Li_' + str(station_index[i]) + '.csv'
        Rrs_path = dir_write + str(station_index[i]) + '_' + str(stations[i]) + '/' + 'reflectance_' + str(station_index[i]) + '.csv'
        exLwn_path = dir_write + str(station_index[i]) + '_' + str(stations[i]) + '/' + 'nLw_' + str(station_index[i]) + '.csv'
        
        write_spectra(Es_OLCI[1], Es_info, Es_path)  # Es_OLCI[1] etc, is pandas dataframe format
        write_spectra(Lt_OLCI[1], Lt_info, Lt_path)
        write_spectra(Li_OLCI[1], Li_info, Li_path)
        write_spectra(Rrs_OLCI[1], Rrs_info, Rrs_path)
        write_spectra(exLwn_OLCI[1], exLwn_info, exLwn_path)
        

        # perform averaging/variability for each spectra at each station -outputs dataframe that is appended to station summary in write function
        Es_av = station_averages_dataframe(Es_OLCI[0],'Es', station_index[i], time_L1, windspeed_L1, phi_L1, srf_bands)
        Lt_av = station_averages_dataframe(Lt_OLCI[0],'Lt', station_index[i], time_L1, windspeed_L1, phi_L1, srf_bands)
        Li_av = station_averages_dataframe(Li_OLCI[0],'Li', station_index[i], time_L1, windspeed_L1, phi_L1, srf_bands)
        Rrs_av = station_averages_dataframe(Rrs_OLCI[0],'reflectance', station_index[i], time_L1, windspeed_L1, phi_L1, srf_bands) 
        exLwn_av = station_averages_dataframe(exLwn_OLCI[0],'nLw', station_index[i], time_L1, windspeed_L1, phi_L1, srf_bands)
        
        Es_av_summary.append(Es_av)
        Lt_av_summary.append(Lt_av)
        Li_av_summary.append(Li_av)
        Rrs_av_summary.append(Rrs_av)
        exLwn_av_summary.append(exLwn_av)
        
        summary_path = dir_write + 'FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_V6_rhoM99_nLw_localChl_NIRcorrection.csv'
        dict_list = append_to_summary(dict_list) # for conditions summary
        
    
    # overall summary - used for station selection    
    write_station_summary(summary_path, time, Es_av_summary, Lt_av_summary, Li_av_summary, Rrs_av_summary, exLwn_av_summary) 
    
    # overall summary - used for station selection    
    ov_summary_path = dir_write +  '/' + 'FICE_conditions_summary.csv'
    df_overall = pd.DataFrame.from_dict(dict_list)  
    df_overall.to_csv(ov_summary_path, na_rep ='NaN', index = False)
        
    # conditions_summaryplots(df_overall)


  