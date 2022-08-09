#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:01:49 2022

@author: tjor

Python implementation of HSAS radiometic post-processing & data formatting for 
FICE 2022

Script includes:

(i) Loading of L1 (Es, Li, Lt) and L2 (Rrs) output .mat files and end-end 
uncertianty estimates

(ii) Spectral downsampling of Li, Lt, Es, Rrs using OLCI spectral response function

(iii) Writing output spectra (downsampled to OLCI bands) to csv file (for each station)
- Includes Rrs, Es, Li, Lt, & metadata for each timestamp 
- Metadata includes: rho, wind-speed, `cloudiness index: pi*Li(400)/Es(400), 
mask for glint (Lt-based) QC, mask for tilt QC, mask for,

(iv) Station-averaging (mean, median, STD) of radiometric and wind data

(v) Writing station summary to .csv file based on FICE_2022_AAOT template
- Includes station avergaes for Rrs, Es, Lt, Lt, Wind & core metadata required
-
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
    srf_bands =  srf_data['nominal_centre_wavelength'][0:16]  # band centres where Hsas has data
    # band_width = srf_data['bandwidth_fwhm'][0:16] - not needed
    srf_wv = np.array(srf_data['mean_spectral_response_function_wavelength'])
    srf = np.array(srf_data['mean_spectral_response_function'])

    return srf, srf_wv, srf_bands


def unpack_L1(fn_L1):
    ''' unpacks relevant L1 data in np array format: L1['L1'].__dict__.keys()
    can be used to show fields - needs to be applied at each struct level'''
    
    L1 = io.loadmat(fn_L1[0], squeeze_me = True, struct_as_record = False) 
    
    time = np.array(numeric_to_UTC_time(L1['L1'].gps.time)) 
    windspeed = L1['L1'].appwind_spd
 
    Es= L1['L1'].instr.Es.data # hyperspectral data
    Lt = L1['L1'].instr.Lt.data
    Li = L1['L1'].instr.Li.data
    
    wv = L1['L1'].wv 

    return time, windspeed, Es, Lt, Li, wv


def unpack_L2(fn_L2, time_bins , wv_bins):
    ''' unpacks relevant L2 data in np array format. Adjustment is made so that 
    size of output matches L1 data. Timestamps that not pass quality control (
    (Lt glint filter) are filled with NaNs'''
  
    L2 = io.loadmat(fn_L2[0], squeeze_me = True, struct_as_record = False) 
    
    # load files in L2 format
    Rrs_L2 = L2['L2'].Rrs.data 
    rho_L2 = L2['L2'].rho 
    exLwn_L2 = L2['L2'].exLwn.data
    qc_ind = L2['L2'].glint_filter.indices_passed_qc  # indices that were kept from L1 data

    # convert to length of L1 data, with nan-padding added where no data
    qc_mask = np.zeros(time_bins) # 1 == data, 0 === no data
    Rrs = np.nan*np.ones([time_bins, wv_bins])
    exLwn = np.nan*np.ones([time_bins, wv_bins])
    rho = np.nan*np.ones(time_bins)
    print(len(qc_ind)) # note - qc ind does not currently match length of Rrs - issue here
    for i in range(len(qc_ind)): # fill L1-size stuctures
        print(i)
        qc_mask[int(qc_ind[i])] = 1 
        Rrs[int(qc_ind[i]),:] = Rrs_L2[i,:] 
        exLwn[int(qc_ind[i]),:] = exLwn_L2[i,:] 
        rho[int(qc_ind[i])] = rho_L2[i]
        
    return qc_mask, Rrs, exLwn, rho


def numeric_to_UTC_time(time_num):
     '''Function to convert numeric time (octave output for Hsas spectra) to UTC 
     - based on Junfangs' uncertainty code'''
   
     time_UTC = []
     for i in range(len(time_num)):
         time_UTC_i = (datetime.datetime.fromordinal(int(time_num[i])) + 
                        datetime.timedelta(seconds = 86400*(time_num[i] - float(int(time_num[i])))) - datetime.timedelta(days=365))
         
         time_UTC.append(time_UTC_i)
     
     return time_UTC


def hyperspec_to_OLCIbands(S, time, wv, srf, srf_bands, srf_wv):
    ''' Function to downsample hyperspectral data to OLCI bands: input spectra 
   in dimenions of time * wavelength'''
   
    D = np.nan*np.zeros([len(time),len(srf_bands)]) # matrix format for down-sampled spectra - time * spectral band

    # Loop for spectral down-sampling (outputs np array D)
    for j in range(len(S)): # loop over timestamps
        wv_j = wv[np.isnan(S[j,:]) == 0] # remove nan padding (rq. for interpolation functions)
        S_j = S[j, np.isnan(S[j,:]) == 0]    
        for k in range(len(srf_bands)): # loop over spectral bands
            interp_funct = interp1d(wv_j, S_j, kind = 'cubic') # interpolate to same wv interval as OLCI
            S_j_interp = interp_funct(srf_wv[k]) # oversample on same wv range as OLCI SRF 
            S_j_k = np.sum(S_j_interp*srf[k])/np.sum(srf[k]) # convolution of kth band for jth timestamp with OLCI SRF
            D[j,k] =  S_j_k # 
    
    # Loop for pandas dataframe format
    D_df = pd.DataFrame(index=time) # pandas data frame format for down-sampled spectra
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


def station_averages_dataframe(D, spec_id, station_i, time, windspeed, srf_bands):
    '''computes averages in data frame format including station meta data '''
          
    summary_df = pd.DataFrame(index=None) 
    summary_df['spectrum'] = [spec_id + '_median', spec_id + '_standard_deviation', spec_id + '_mean', spec_id + '_uncertainty'] 
    summary_df['station'] = [station_i, station_i, station_i, station_i]
    summary_df['water body'] = ['Adriatic Sea', 'Adriatic Sea', 'Adriatic Sea', 'Adriatic Sea']
    summary_df['station start time'] = [str(time[0])[0:19], str(time[0])[0:19], str(time[0])[0:19], str(time[0])[0:19]]
    summary_df['wind speed [m/s] (median, std, mean)'] = [str(np.nanmedian(windspeed)), str(np.nanstd(windspeed)), str(np.nanmean(windspeed)), '']
    summary_df['sechi depth'] = ['', '', '', '']  
    
    # averages and std and uncertainty  of downsampled spectra in each band (np.array format)
    D_summary = np.zeros([4,len(srf_bands)])
    D_summary[0,:] = np.nanmedian(D,axis=0)
    D_summary[1,:] = np.nanstd(D,axis=0)
    D_summary[2,:] = np.nanmean(D,axis=0)
    D_summary[3,:] = np.zeros(len(srf_bands))  # temporary fix for uncertainty

    # fill data frame with np.array info
    for k in range(len(srf_bands)):
        summary_df[str(srf_bands[k])] = D_summary[:,k] 
  
    return summary_df


def write_spectra(D_df,sensor_info, fpath):
    ''' writes downsampled spectra to csvfile'''
    
    # create headers using sensor info dictionary
    with open(fpath, 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        for i in range(len(sensor_info)):
            row_i = str('# '  + list(sensor_info.keys())[i] + ': ' + str(list(sensor_info.values())[i]))
            print(row_i)
            writer.writerow([''.join(row_i)])
    
    # append spectra dataframe to csvfile: dim = time * OLCI bands
    D_df.to_csv(fpath, mode = 'a')
    
    return 

       
def write_station_summary(summary_path, time, Es_av, Lt_av, Li_av):
    ''' Function to write station summary in csv format'''
    
    # fill up pre-amble
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
        row = str('#time_period: [{"start:"' + str(time[0])[0:19] + ',"end:"' + str(time[-1])[0:19] + '}]')
        writer.writerow([''.join(row)])
        row = str('#reflectance_calculation equation: TO ADD')
        writer.writerow([''.join(row)])
        row = str('units: {"reflectance":["Wavelength [nm]","Reflectance [ / sr"],"es":["Wavelength [nm]","Irradiance [mW / (cm^2 * 10^-6 m)]"],"lt":["Wavelength [nm]","Radiance [mW / (sr * cm^2 * 10^-6 m)]"],"li":["Wavelength [nm]","Radiance [mW / (sr cm^2 * 10^-6 m)]"]}')
        writer.writerow([''.join(row)])
    
    # append station averages to csv file - refletance and nLW to add latter
    Es_av.to_csv(summary_path, index = False, mode = 'a')
    Lt_av.to_csv(summary_path, index = False, mode = 'a')
    Li_av.to_csv(summary_path, index = False, mode = 'a')
    
    return

    
if __name__ == '__main__':


    # subdirectories to read data
    dir_main = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/'
    dir_L0 = dir_main + 'Processed_TPcorrection/L0/'
    dir_L1 = dir_main + 'Processed_TPcorrection/L1/'
    dir_L2 = dir_main + 'Processed_TPcorrection/L2/'
    dir_unc = dir_main + 'TO_ADD' 
    dir_srf = dir_main + '/Source/HSAS/HSAS_Processing/FICE_Formatting' # cuntains OLCI SRF and template
    
    # subdirectory to write data
    dir_write = dir_main +'Processed_TPcorrection/FICE_format/'

    # load OLCI srf
    srf, srf_wv, srf_bands = unpack_srf(dir_srf)

    # initialize station sub directories
    stations =  os.listdir(dir_L1) # each subdirectory is a station
    # for i in range(len(stations)):
       # os.mkdir(dir_write + str(stations[i]))

    for i in range(3): # process each station in sequence
    
        # access mat filenames (hsas data structures) of ith station 
        # fn_L0_i = glob.glob(dir_L0 + stations[i]  + '/*mat*') 
        fn_L1 = glob.glob(dir_L1 + stations[i]  + '/*mat*') # applies to data from ith station
        fn_L2 = glob.glob(dir_L2 + stations[i]  + '/*mat*')
        
        # unpack variables from L1 and L2 data data stuctures in np array format
        time, windspeed, Es, Lt, Li, wv = unpack_L1(fn_L1) # applies to data from ith station 
        # qc_mask, Rrs, exLwn, rho = unpack_L2(fn_L2, len(time) , len(wv)) # applies to data from ith station
        
        # spectral downsampling to  OLCI: element 0 is np array, 1 is dataframe formate
        Es_OLCI = hyperspec_to_OLCIbands(Es, time, wv, srf, srf_bands, srf_wv)  
        Lt_OLCI = hyperspec_to_OLCIbands(Lt, time, wv, srf, srf_bands, srf_wv)
        Li_OLCI = hyperspec_to_OLCIbands(Li, time, wv, srf, srf_bands, srf_wv)
        # multi_hyper_plot(Es, wv, Es_OLCI, srf_bands, 'Es [mW cm$^{2}$ nm$^{1}$')
        
        # collate sensor info and write downsampled spectra to csv files
        Es_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '2027A', 'Integration time': 'TO ADD', 'Raw wavelength scale': 'TO ADD'} 
        Lt_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '464', 'Integration time': 'TO ADD', 'Raw wavelength scale': 'TO ADD'} 
        Li_info = {'Make': 'Seabird', 'Model': 'HypserSAS', 'Serial number': '2054A', 'Integration time': 'TO ADD', 'Raw wavelength scale': 'TO ADD'} 
        
        Es_path = dir_write +  stations[i] + '/' + 'Es_' + stations[i] + '.csv'
        Lt_path = dir_write +  stations[i] + '/' + 'Lt_' + stations[i] + '.csv'
        Li_path = dir_write +  stations[i] + '/' + 'Li_' + stations[i] + '.csv'
        
        write_spectra(Es_OLCI[1], Es_info, Es_path)  # Es_OLCI[1] is pandas dataframe format
        write_spectra(Lt_OLCI[1], Lt_info, Lt_path)
        write_spectra(Li_OLCI[1], Li_info, Li_path)
        
        #Es_OLCI_mean, Es_OLCI_med, Es_OLCI_std = station_averages(Es_OLCI[0]) #  Es_OLCI[0] is np array format
        #Lt_OLCI_mean, Lt_OLCI_med, Lt_OLCI_std = station_averages(Lt_OLCI[0])
        #Li_OLCI_mean, Li_OLCI_med, Li_OLCI_std = station_averages(Li_OLCI[0])
     
        Es_av = station_averages_dataframe(Es_OLCI[0],'Es', stations[i], time, windspeed, srf_bands)
        Lt_av = station_averages_dataframe(Lt_OLCI[0],'Lt', stations[i], time, windspeed, srf_bands)
        Li_av = station_averages_dataframe(Li_OLCI[0],'Li', stations[i], time, windspeed, srf_bands)
        
        summary_path = dir_write +  stations[i] + '/' + 'Summary_' + stations[i] + '.csv'
        write_station_summary(summary_path, time, Es_av, Lt_av, Li_av)
       
        
        

 