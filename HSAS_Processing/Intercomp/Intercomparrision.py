#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:36:24 2022

@author: tjor

This script performs the system inter


# notes 
    
# 1. Baselines

 Ed, Lt, Ls, Rrs: use 4-way mean of NASA, PML, Tartu, HEREON
 nLw: use SEAPRISM-Aeronet as independnet reference(Zibordi et al. 2018)

# 2. QC levels
 All data
 Zibordi's Aeronet QC 
 Cloudiness filtering (subset passed)
 Cloudiness filtering (subset failed) ?? not needed ??
 
# 3. Methodology follows Tilstone 2020.

    Match-up stats - bands 443, 560 and 665 nm
        N - no. of stations used for 
        RPD - relative % difference
        RMSD - root mean square diviation
        
    Scatter plots and spectral depedence of residuals as a box plot 
    (show all bands simulatneously)

#

"""

from csv import reader
import numpy as np
import pandas as pd
import datetime 
import os

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches

import dataframe_image as dfi

######################
# READ FUNCTIONS
#######################

# PML read functions
def read_PML_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_PML_spec(path, 'Es_mean') 
    Lsky = read_PML_spec(path, 'Li_mean') 
    Lt = read_PML_spec(path, 'Lt_mean')
    Rrs = read_PML_spec(path, 'reflectance_mean') 
    Rrs_std = read_PML_spec(path, 'reflectance_standard_deviation') 
    nLw = read_PML_spec(path, 'nLw_mean') 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw

def read_PML_spec(path, S): 
    '''Returns dataframe for spectral type S, OLCI bands as columns, indexed by station
    ''''''PML has slightly different read function than other systems as we require all timestamps 
    (used as reference time)'''

    
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)
        
        station = np.arange(0,78)
        time_start = np.empty(78, dtype=object)
        windspeed = np.nan*np.ones(78)
        spec_data = np.nan*np.ones([78,19]) 
        azimuth = np.nan*np.ones(78)
    
        # fill up data matrix
        i = 0
        for row in csv_reader:
            if i == 70:
                station[i] = 70
                time_start[i] = '2022-07-20 14:00:00' # missing staion- TS added (from Tartu)
                windspeed[i] = 2.7
                azimuth[i] =  90
                i=i+1
            if(row[0]) == S: 
                    station[i] = row[1]
                    if station[i] != 2: # cleaning for duplicate
                        timestamp_i =  datetime.datetime.strptime(row[3], '%Y-%m-%d %H:%M:%S')
                        time_start[i] = str(timestamp_i)
                        windspeed[i] = row[4]
                        azimuth[i] = row[5]
                        spec_data[i,:] = row[6:25]
                    else:
                        time_start[i]= '2022-07-14 09:00:00' # duplicate station (from Tartu)
                        windspeed[i] = 1.2
                        azimuth[i] =  135
                    i = i + 1

        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        for i in range(len(bands)):
                df[str(bands[i])] = spec_data[:,i] 
       
    return df
    

# NASA read functions 
def read_NASA_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_NASA_spec(path, 'ed_average') 
    Lsky = read_NASA_spec(path, 'lsky_average') 
    Lt = read_NASA_spec(path, 'lu_average')
    Rrs = read_NASA_spec(path, 'reflectance_average') 
    Rrs_std = read_NASA_spec(path, 'reflectance_uncertainty')#  - unc not std?
    nLw = read_NASA_spec(path, 'nlw_average') 

    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw

def read_NASA_spec(path, S): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        data = np.nan*np.ones([78,19]) 
        time_start = np.empty(78, dtype=object)
        windspeed = np.nan*np.ones(78)
        azimuth = np.nan*np.ones(78)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                i = int(row[1]) 
                timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                time_start[i] = str(timestamp_i)
                windspeed[i] = row[4]
                azimuth[i] = row[5]
                data[i, :] = row[6:25]
             
         
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        
        for i in range(len(bands)):
                if S =='ed_average' or S == 'lsky_average' or S == 'lu_average' or S == 'nlw_average':
                    df[str(bands[i])] = 10*data[:,i] 
                else:
                    df[str(bands[i])] = data[:,i] 
   
    return df


# Tartu read functions 
def read_TARTU_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_TARTU_spec(path, 'ramses3_above_ed_average') 
    Lsky = read_TARTU_spec(path, 'ramses3_above_lsky_average')
    Lt = read_TARTU_spec(path, 'ramses3_above_lu_average')
    Rrs = read_TARTU_spec(path, 'ramses3_above_reflectance_average') 
    Rrs_std = read_TARTU_spec(path, 'ramses3_above_reflectance_standardDeviation')
    nLw = read_TARTU_spec(path, 'ramses3_above_nlw_average') # yet to arrive
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw

def read_TARTU_spec(path, S): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    'Timestamp is read in np64 format'
    
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)
        
        station = np.arange(0,78)
        time_start = np.empty(78,  dtype=object)
        windspeed = np.nan*np.ones(78)
        spec_data = np.nan*np.ones([78,19]) 
 
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                i = int(row[3]) 
                timestamp_i = datetime.datetime.strptime(row[5][0:19], '%Y-%m-%dT%H:%M:%S')
                time_start[i] = str(timestamp_i)
                windspeed[i] = row[6]
                spec_data[i,:] = row[10:29]
   

        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        for i in range(len(bands)):
            if S ==  'ramses3_above_reflectance_average':
                df[str(bands[i])] = spec_data[:,i]/np.pi
            elif  S ==  'ramses3_above_reflectance_standardDeviation':
                df[str(bands[i])] = spec_data[:,i]/np.sqrt(np.pi)
            else:
                df[str(bands[i])] = spec_data[:,i]
        
    return df


# RBINS read functions 
def read_RBINS_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_RBINS_spec(path, 'ed_average') 
    Lsky = read_RBINS_spec(path, 'lsky_average')
    Lt = read_RBINS_spec(path, 'lu_average')
    Rrs = read_RBINS_spec(path, 'reflectance_average') 
    Rrs_std = read_RBINS_spec(path, 'reflectance_standardDeviation') 
    nLw = read_RBINS_spec(path, 'nlw_average') 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw


def read_RBINS_spec(path, S): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        spec_data = np.nan*np.ones([78,19]) 
        time_start = np.nan*np.ones(78,dtype=object)
        windspeed = np.nan*np.ones(78)
        azimuth = np.nan*np.ones(78)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                if str(row[1][0]) != 'P': #
                    i = int(row[1][0:2])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    time_start[i] = str(timestamp_i)
                    windspeed[i] = row[4]
                    azimuth[i] = 135
                    spec_data[i,:] = row[6:25]

            
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        for i in range(len(bands)):
            if S == 'reflectance_average':
                df[str(bands[i])] = spec_data[:,i]/np.pi
            elif  S ==  'ramses3_above_reflectance_standardDeviation':
                df[str(bands[i])] = spec_data[:,i]/np.sqrt(np.pi)
            else:
                df[str(bands[i])] = spec_data[:,i]
  
    return df


# HEREON read functions 
def read_HEREON_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_HEREON_spec(path, 'ed_average') 
    Lsky = read_HEREON_spec(path, 'lsky_average')
    Lt = read_HEREON_spec(path, 'lt_average')
    Rrs = read_HEREON_spec(path, 'reflectance_average') 
    Rrs_std = read_HEREON_spec(path, 'reflectance_standardDeviation') 
    nLw = read_HEREON_spec(path, 'lw_average') 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw


def read_HEREON_spec(path, S): 
    ' read function for HEREON FILE FORMAT: First Line: Wavelength; 1) Rrs_median, 2) Rrs_std, 3) Rrs_average, 4) Ed_median, 5) Ed_std, 6) Ed_average, 7) Lu_median, 8) Lu_std, 9) Lu_average, 10) Lsky_median, 11) Lsky_std, 12) Lsky_average, 13) Lw_median, 14) Lw_std, 15) Lw_average'
   
    files = os.listdir(path)
    data = np.nan*np.ones([78,19]) 
    station = np.arange(0,78)
       
    if   S == 'ed_average':
        row_index = 6
    elif S == 'lt_average':
        row_index = 9
    elif S == 'lsky_average':
        row_index = 12
    elif S == 'reflectance_average':
        row_index = 3
    elif S == 'reflectance_standardDeviation':
            row_index = 2
    elif S == 'lw_average':
        row_index = 15

    for i in range(len(files)):
        station[i] = i
        text_i = np.loadtxt(path_HEREON + '/' + files[i], skiprows = 1)
        data[i,:] = text_i[row_index,1:-1] 
        
        df = pd.DataFrame(index = station) 
        for i in range(len(bands)):
                df[str(bands[i])] = data[:,i] 
    
    return df

# NOAA read functions 
def read_NOAA_data(path, bands, Ed_PML):
    'Wrappers to read spectra in dataframe format - Ed_PML used for reference time'
    
    Ed = read_NOAA_spec(path, 'Es_average') 
    Lsky = read_NOAA_spec(path, 'Lsky_average') 
    Lt =  read_NOAA_spec(path, 'Lt_average') 
    Rrs =  read_NOAA_spec(path, 'Rrs_average') 
    nLw = read_NOAA_spec(path, 'nlw_average')

    Ed = time_match_NOAA(Ed, Ed_PML)
    nLw = time_match_NOAA(nLw, Ed_PML)

    return Ed, Lsky, Lt, Rrs, nLw


def read_NOAA_spec(path, S): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    
    with open(path, 'r') as read_obj:
       
        # fill up data matrix  - case ed
        if  S == 'Es_average':
                
            csv_reader = reader(read_obj)  
            data = np.nan*np.ones([78,16]) 
            time_start = np.empty(78, dtype=object)
            # station = np.arange(0,78)
            i = 0
            for row in csv_reader:
                if(row[0]) == S: 
                    # print(row[6:22])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    time_start[i] = str(timestamp_i)
                    data[i,:] = row[6:22]
                    i = i +1
                
            # convert to df format - 
            df = pd.DataFrame() 
            df['time_start'] = time_start
            
            for i in range(len(bands) - 3): # considers bands < 800 nm              
                df[str(bands[i])] = 10*data[:,i] 
           
       # fill up data matrix - case nLw or Lw
        elif S == 'nlw_average' or S == 'lw_average':
                        
            csv_reader = reader(read_obj)  
            data = np.nan*np.ones([78,8]) 
            time_start = np.empty(78, dtype=object)
            # station = np.arange(0,78)
            i = 0
            for row in csv_reader:
                if(row[0]) == S: 
                    print(row[6:14])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    time_start[i] = str(timestamp_i)
                    data[i, :] = row[6:14]
                    i = i +1
                
            # convert to df format - 
            df = pd.DataFrame() 
            df['time_start'] = time_start
            
            for i in range(8): # considers bands < 800 nm              
                df[str(bands[i])] = 10*data[:,i] 
                
        else:
            df = pd.DataFrame() 
            
    return df

def time_match_NOAA(df_NOAA, Ed_PML):
    'Sub routine to find station index for NOAA, and then reformat index of dataframe - assumes 10 min buffer'
    
    time_start_PML = [datetime.datetime.strptime(Ed_PML['time_start'][i],'%Y-%m-%d %H:%M:%S')  for i in range(len(Ed_PML))] 
    time_start_NOAA = [datetime.datetime.strptime(df_NOAA['time_start'][i],'%Y-%m-%d %H:%M:%S')  for i in range(72)]        # NOAA just has 72 timestamps     
    
    time_start_NOAA_matching = np.nan*np.ones(78,dtype = object)
    spec_data_matching = np.nan*np.ones([len(df_NOAA),len(df_NOAA.columns)]) 
    tol = 10*60       
    for i in range(len(time_start_PML)):
         nearest_time, nearest_index = nearest(time_start_NOAA, time_start_PML[i])
         delta_t = abs(time_start_PML[i] - nearest_time) 
         print(delta_t)
         print(nearest_index)
         if delta_t.total_seconds() < tol:
             time_start_NOAA_matching[i] = str(time_start_NOAA[nearest_index]) 
             for j in range(len(df_NOAA.columns)-1):
                 spec_data_matching[i,j] = df_NOAA[str(bands[j])][nearest_index]
      
    # convert to df format - columns are station, time, windspeed, spectra in each band
    df = pd.DataFrame() 
    df['time_start'] =  time_start_NOAA_matching
    for j in range(len(df_NOAA.columns)-1): # considers bands < 800 nm              
        df[str(bands[j])] =  spec_data_matching[:,j] 
           
    return df

# CNR read functions 
def read_CNR_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_CNR_spec(path, 'ed_average') 
    Lsky = read_CNR_spec(path, 'lsky_average') 
    Lt = read_CNR_spec(path, 'lu_average')
    Rrs = read_CNR_spec(path, 'reflectance_average') 
    Rrs_std = read_CNR_spec(path, 'reflectance_stdev')
    nlw = read_CNR_spec(path, 'nlw') # - unc not std?

    return Ed, Lsky, Lt, Rrs, Rrs_std, nlw

def read_CNR_spec(path, S): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        data = np.nan*np.ones([78,19]) 
        time_start = np.empty(78, dtype=object)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                station = row[1]
                if station != 'None':
                    i = int(row[1]) 
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%MZ')
                    time_start[i] = str(timestamp_i)
                    data[i, :] = row[6:25]
             
         
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame() 
        df['time_start'] = time_start
        if S == 'reflectance_average' or S == 'reflectance_stdev':
            for i in range(len(bands)):
                  df[str(bands[i])] = data[:,i]/np.pi 
        else:
            for i in range(len(bands)):
                  df[str(bands[i])] = data[:,i] 
                  
             
    return df

####################
# BASELINES
#####################

def baseline_average(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON):
    'Computes reference baseline (_R) based on weighthed mean of PML, NASA, HEREON, TARTU.'                             
    
    if spec_type == 'Ed':
        # first, average seabird reference systems
        df_seabird = ( 
            # combine dataframes into a single dataframe
            pd.concat([df_PML, df_NASA])
            # replace 0 values with nan to exclude them from mean calculation
            .replace(0, np.nan)
            .reset_index()
            # group by the row within the original dataframe
            .groupby("index")
            # calculate the mean
            .mean()
        )
        
    elif spec_type == 'Lsky' or  spec_type == 'Lt' or spec_type == 'Rrs':
        # first, average seabird reference systems
        df_seabird = df_PML
       
    # second, average trios reference systems
    df_trios = (
        # combine dataframes into a single dataframe
        pd.concat([df_TARTU, df_HEREON])
        # replace 0 values with nan to exclude them from mean calculation
        .replace(0, np.nan)
        .reset_index()
        # group by the row within the original dataframe
        .groupby("index")
        # calculate the mean
        .mean()
    )
    
    # average trios and seabird
    df_R = (
        # combine dataframes into a single dataframe
        pd.concat([df_seabird, df_trios])
        # replace 0 values with nan to exclude them from mean calculation
        .replace(0, np.nan)
        .reset_index()
        # group by the row within the original dataframe
        .groupby("index")
        # calculate the mean
        .mean()
    )
    
    # used to calculate number used in mean - require need 3 or more systems for mean to be defined
    df_N = pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
    N_mean = np.zeros(78)
    for i in range(len(df_R)):
        no_of_records = 0
        for j in range(4):
            if np.isnan(df_N['400'][i].iloc[j]) == False:
                no_of_records =  no_of_records + 1
            N_mean[i] = no_of_records
            
    # PML used for timestamps and windspeeds (PML is complete dataset  so copied to reference)    
    df_R['time_start']  = df_PML['time_start'] 
    df_R['windspeed']  = df_PML['windspeed']     
    df_R['N_mean']  = N_mean    
    
    return df_R


def baseline_average_V2(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON):
      'Computes reference baselines - assummes even weighting of each submitted team'''
    
      df_R = ( 
            # combine dataframes into a single dataframe
            pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
            # replace 0 values with nan to exclude them from mean calculation
            .replace(0, np.nan)
            .reset_index()
            # group by the row within the original dataframe
            .groupby("index")
            # calculate the mean
            .mean()
        )
       
       # used to calculate number used in mean - require need 3 or more systems for mean to be defined
      df_N =  pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
      N_mean = np.zeros(78)
      for i in range(len(df_R)):
           no_of_records = 0
           for j in range(4):
               if np.isnan(df_N['400'][i].iloc[j]) == False:
                   no_of_records =  no_of_records + 1
               N_mean[i] = no_of_records
               
       # PML used for timestamps and windspeeds (PML is complete dataset  so copied to reference)    
      df_R['time_start']  = df_PML['time_start'] 
      df_R['windspeed']  = df_PML['windspeed']     
      df_R['N_mean']  = N_mean    
         
      return df_R


def read_Aeronet_nLw(path_NLW, Ed_R, bands):
    ''' function to read nLw 1.5 data - assumes similar format to baseline _R dataframes
    Based on nearest neighbour interpolation within a +/- 10 minute tolerance to be counted.
    Note : 667 nm band is currently interpreted as 665 nm OLCI band (may need interpolating later on)
    Units are mw/(cm^2 sr micro-m) - multiply by 10 to get  mw/(m^2 sr nm)'''
    
    # extract nLw data in data frame formate
    files = os.listdir(path_NLW)
    data = pd.DataFrame() 
    for i in range(len(files)): 
        data_i = pd.read_csv(path_NLW + '/' + files[i], skiprows = 0)
        data.append(data_i)
        data = pd.concat([data, data_i]) #
    data_AOC = data.reset_index(drop=True)
    
    # convert timestamps to datetime formats - Ed_R is used for station timestamps
    date_nLw = [datetime.datetime.strptime(str(data_AOC['Date(dd-mm-yyyy)'][i]),'%d:%m:%Y') for i in range(len(data_AOC))] # combine dates and times ino single timestamp
    time_nLw = [datetime.datetime.strptime(str(data_AOC['Time(hh:mm:ss)'][i]),'%H:%M:%S') for i in range(len(data_AOC))]
    timestamp_nLw = [datetime.datetime.combine(date_nLw[i].date(), time_nLw[i].time()) for i in range(len(data_AOC))] 
    timestamp_station = [datetime.datetime.strptime(Ed_R['time_start'][i],'%Y-%m-%d %H:%M:%S') + datetime.timedelta(0,2,30) for i in range(len(Ed_R))]       
         
    # fill up nLw data matrix based on time stamp matchung (similar approach to QC) - 
    keys_nLw = ['Lwn_f/Q[412nm]','Lwn_f/Q[443nm]','Lwn_f/Q[490nm]','Lwn_f/Q[510nm]','Lwn_f/Q[560nm]','Lwn_f/Q[620nm]','Lwn_f/Q[667nm]']
    nLw_data = np.nan*np.zeros([78, len(keys_nLw)])
    station = np.arange(0,78,1)
    tol = 10*60 # as with QC flags +/- 10 minutes is default tolerance (should approximately match QC)
    for i in range(len(timestamp_station)):# Loop over FICE timestamps to find good AOC flag with +/- tolerance
        nearest_time, nearest_index = nearest(timestamp_nLw, timestamp_station[i])
        delta_t = abs(timestamp_station[i] - nearest_time) # time difference
        if delta_t.total_seconds() < tol:
            for j in range(len(keys_nLw)):
                nLw_data[i,j] = 10*data_AOC[keys_nLw[j]][i] # factor 10 for unit conversion
                
    # fill up data frame with similar format as baseline references
    df_R = pd.DataFrame(index = station) 
    for j in range(len(keys_nLw)):
            df_R[str(bands[j+1])] = nLw_data[:,j]
    df_R['time_start']  = Ed_R['time_start'] 
    
    return df_R
    

# QC functions
def nearest(items, pivot): 
    'function to locate nearest values and indicies in list x relative to a pivot (used to locate nearest timestamps)'
    
    nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value    
    nearest_index= items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index
   
    return nearest_value, nearest_index


def QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML):
    
    ''' Function to derive QC mask considering: (i) Aeronet QC, (ii) Filtering based
    on cloudiness and CV_rrs (averaged on OLCI bands 1-6). 1==good data, 0==bad data
   
    - Method (i) icludes centering data timestamp to middle of 5 min recording interval.
    The threshold rank_AOC = AOC['Rank'] # rank 0.6 or 1 is a QC pass
   
    - Ed_R, Lsky_R, Rrs_R, Rrs_std_R  (4-way reference values) are used in filtering step 
    for method (ii).
    
    - tol is the maximum time (+/-) that a good Aeronet OC timestamp can be from a FICE 
    data record for the fICE
    

    '''   
    
    # (i) QC using Aeronet QC flags from Zibordiplt.figure(figsize=(12,12))
    
    # Read AOC data and convert FICE timestamp to dt format
    AOC_data = pd.read_csv(path_QC, delim_whitespace=True, skiprows=[0])   # QC mask from aeronet
    timestamp_AOC = [datetime.datetime.strptime(AOC_data['Lev_1.5_record_ID'][i][5:24],'%d:%m:%Y-%H:%M:%S') for i in range(len(AOC_data))]  
    rank_AOC = AOC_data['Rank'] 
    timestamp_data = [datetime.datetime.strptime(Ed_R['time_start'][i],'%Y-%m-%d %H:%M:%S') + datetime.timedelta(0,2,30) for i in range(len(Ed_R))] 
    
    # Loop over FICE timestamps to find good AOC flag with +/- tolerance
    tol = 10*60 # +/- tolerance in seconds - +/1 10 minutes used as default
    QC_AOC = np.zeros(78) # AOC_QC mask: zero == failed QC, one = passed QC
    QC_AOC_3 = np.zeros(78)  # AOC_QC with at least 3 systems present
    for i in range(len(timestamp_data)):
      nearest_time, nearest_index = nearest(timestamp_AOC, timestamp_data[i]) # finds nearest AOC timestamp and index to FICE record
      delta_t = abs(timestamp_data[i] - nearest_time) # time difference
      if delta_t.total_seconds() < tol and rank_AOC[nearest_index] >= 0.6: # test if time difference and rank criteria are satisfied
          QC_AOC[i] = 1 # 
      if delta_t.total_seconds() < tol and rank_AOC[nearest_index] >=0.6 and Ed_R['N_mean'][i] > 2:
          QC_AOC_3[i] = 1 # 
    # (ii) Cloudiness index QC  
    
    CI = np.pi*Lsky_PML['400']/Ed_PML['400'] # cloudiness index from Simis 2013. Less than  0.4 is clear or light haze
    CV_Rrs_specav = (1/6)*100*(Rrs_std_PML['400']/Rrs_PML['400'] + Rrs_std_PML['412.5']/Rrs_PML['412.5']+ Rrs_std_PML['442.5']/Rrs_PML['442.5']+ Rrs_std_PML['490']/Rrs_PML['490']+ Rrs_std_PML['510']/Rrs_PML['510'] + Rrs_std_PML['560']/Rrs_PML['560'] + Rrs_std_PML['560']/Rrs_PML['560']) # baverage on bands 1-6
    
    #
    plt.figure(figsize=(10, 8))
    plt.rc('font', size=16)   
    plt.title('Dependence of Rrs variability on cloudiness (station average values)')
     # colors = cm.jet(np.linspace(0, 1, len(df_overall)))
    #for i in range(len(df_overall)):
    plt.scatter(CV_Rrs_specav,CI,color='blue')
    #  plt.legend(labels = df_overall['station'], ncol=2,fontsize=18)
    plt.xlabel('CV[$R_{rs}$]: mean for OLCI bands 1-6 [%]')
    plt.ylabel('$\pi$$L_{sky}$(400)/$E_{d}$(400): (cloudiness index)')
    plt.ylim(0,1)
    plt.xlim(0,10)

    filename  =  path_output +  '/' + 'Method2QC.png'
    plt.savefig(filename)


    # Loop over FICE timestamps to find good condtions flag
    QC_CI = np.zeros(78) # mask that passes CI QC
    QC_CI_3 = np.zeros(78) # mask that passes CI QC and at least 3 systems present
    for i in range(len(timestamp_data)):
        if CI[i] < 0.4 and CV_Rrs_specav[i] < 2.5:
            QC_CI[i] = 1
        if CI[i] < 0.4 and CV_Rrs_specav[i] < 2.5 and  Ed_R['N_mean'][i] > 2:
            QC_CI_3[i] = 1
            
    
    #  fill up dataframe with QC fields present 
    station = np.array(Ed_R.index)
    df = pd.DataFrame(index = station)  
   # df['time_start'] =  Ed_R['time_start']
    #df['rank_AOC'] = rank_AOC
    #df['QC_AOC'] = QC_AOC
    df['QC_AOC_3'] = QC_AOC_3 # default AOC QC
    #df['Cloudiness_Index'] = CI
    #df['windspeed'] =  Ed_R['windspeed'] 
    #df['CV_Rrs_specAv'] = CV_Rrs_specav
    #df['QC_CI'] = QC_CI # 
    df['QC_CI_3'] = QC_CI_3 # default C.I QC
    
    
    filename  =  path_output +  '/' + 'AAOT_QCmask.csv'
    df.to_csv(filename, na_rep ='NaN')
    
    return df


def plot_dataandQCmasks():
   
    '''plot to show data and QC masks for AAOT: (i) Data no QC, (ii) QC, 
    (iii) Data with AOC QC, (vi) Data with cloudiness index w'''
    # 
    
    
    # 1. plot of all data records
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1: # replace with CNR
                Z[1,i] = 1
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)       
    cmap1 = matplotlib.cm.get_cmap('coolwarm')
    plt.pcolor(Z,edgecolors='k',cmap=cmap1)
    plt.title('AAOT station mask')
    
    plt.xlabel('Station number')
    ylabels = ['NOAA', 'CNR','RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    filename  =  path_output +  '/' + 'stationmask.png'
    plt.savefig(filename)
    
    # 2. plot of QC
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([2, 78])
    Z[0,:] = Q_mask['QC_CI_3']
    Z[1,:] = Q_mask['QC_AOC_3']
  
    # figure for stations
    plt.figure(figsize=(14,6))
    plt.rc('font', size=18)      
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('Quality control mask')
    ylabels = ['2. Env.' + '\n' + '& CV[$R_{rs}$]', '1. AOC:'  +'\n' +  'SEAPRISM']
    plt.xlabel('Station number')
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Pass QC')
    blue_patch = mpatches.Patch(color=cmap1(0), label='Fail QC')
    plt.legend(handles=[red_patch,blue_patch],loc=4)  
    
    filename  =  path_output +  '/' + 'QCmask.png'
    plt.savefig(filename,dpi=300)
    
    # 3. plot of data wwith AOC QC
    
    # creates plotting mask ,Z, where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[1,i] = 1    
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' +  'AOC-SEAPRISM QC and $\geq$ 3 reference systems with data')
    plt.xlabel('Station number')
    ylabels = ['NOAA','CNR', 'RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output +  '/' + 'stationmask_QC1.png'
    plt.savefig(filename)
    
    # 3. plot of data with Cloudiness Index QC

    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[1,i] = 1
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' + 'C.I and CV[$R_{rs}$] QC and $\geq$ 3 reference systems with data')
    
    plt.xlabel('Station number')
    ylabels = ['NOAA','CNR', 'RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output +  '/' + 'stationmask_QC2.png'
    plt.savefig(filename)
    
    return


def scatter_subplot(spec_type, system, plot_index, ylab, xlab, limits, ticks, df_sys, df_R):
    ''' suplot routine for scatter plot'''
   
    colors = cm.rainbow(np.linspace(0,1,10))# color mask to match rrs with time series   
    plt.subplot(2,4,plot_index)    
    X = np.arange(0,2000,1)
    plt.plot(X,X,color='gray')
    plt.title(system)
    
    # fill up
    if spec_type == 'nLw':
        for i in range(1,8,1):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')  
    else :
        for i in range(10):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')
    
    if plot_index == 1: 
        plt.legend(fontsize=10)
    
    # ables
    plt.gca().set_aspect('equal')
    plt.xlim(limits)
    plt.ylim(limits)
    plt.xticks(ticks)
    plt.yticks(ticks)
    
    if plot_index < 4:
        plt.xticks() 

    if spec_type == 'Rrs':    
        plt.xticks(ticks, rotation=45)
        plt.yticks(ticks)
    
    return


def plot_scatter(spec_type,df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, df_CNR, df_NOAA,  Q_mask, Qtype = 'AOC_3'):
    ''' suplot routine for scatter plot'''
  
    # qc filter
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
    df_CNR = df_CNR[Q_mask[Qtype]==1]
    df_NOAA = df_NOAA[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
    
    # scatter plot figure
    fig = plt.figure(figsize=(24,12))
    plt.rc('font', size=18)      
    if spec_type == 'Ed':
        xlab ='Reference: $E_{d}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        ylab = '$E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        limits = [800, 1650]
        ticks = [800, 1000, 1200,1400, 1600]
        plt.suptitle('System inter-comparison for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lsky':
        xlab ='Reference: $L_{sky}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 130]
        ticks = [0, 30, 60, 90, 120]
        plt.suptitle('System inter-comparison for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lt':
        xlab ='Reference: $L_{t}^{r}$(0$^{+})$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{t}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Rrs':
        xlab ='Reference: $R_{rs}^{r}$($\lambda)$ [sr$^{-1}$]'
        ylab = '$R_{rs}$($\lambda)$ [sr$^{-1}$]'
        limits = [0, 0.016]
        ticks = [0, 0.004, 0.008, 0.012, 0.016]
        plt.suptitle('System inter-comparison for remote-sensing reflectance: $R_{rs}$($\lambda$)')
    elif spec_type == 'nLw':
        xlab ='Reference: $L_{wn}^{r}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{nw}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for normalized water-leaving radiance: $L_{wn}$($\lambda$)')
       
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    # subplots
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_PML, df_R)
    
    subtitle = 'HEREON: N = '  + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    scatter_subplot(spec_type,subtitle, index,  ylab, xlab, limits, ticks, df_HEREON, df_R)
    
    subtitle = 'TARTU: N = '  + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_TARTU, df_R)
    
    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_NASA, df_R)

    subtitle = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
    index = 5
    scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_RBINS, df_R)
    
    subtitle = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
    index = 6
    scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_CNR, df_R)
    
    if spec_type =='Ed' or spec_type == 'nLw':
        subtitle = 'NOAA: N = ' + str(np.sum(~np.isnan(df_NOAA['400'])))
        index = 7
        scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_NOAA, df_R)

    plt.tight_layout(pad=1.8)

    
    filename  =  path_output +  '/' + spec_type + '_scattterplot.png'
    plt.savefig(filename)
    
    return


def _resid_subplot(spec_type,system, plot_index, ylab, percent_limits, df_sys, df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,10)) 
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])
    else:
        resid = [] # residual distribution in each band
        for i in range(10):
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i =  resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)

    #
    plt.subplot(2,4,plot_index) 
    plt.title(system)
    plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        bp['boxes'][i].set_facecolor(colors[i])
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands[0:10])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
    
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
     #   plt.xlabel('Wavelength [nm]')
  
    return


def plot_residuals(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, df_CNR, df_NOAA, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
    df_CNR = df_CNR[Q_mask[Qtype]==1]
    df_NOAA = df_NOAA[Q_mask[Qtype]==1]
    
    
    df_R = df_R[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,12))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual [%]'
        plt.suptitle('Percentage residuals for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 6
        percent_limits_2 = 12
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual [%]'
         plt.suptitle('Percentage residuals for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 6
         percent_limits_2 = 12
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual [%]'
         plt.suptitle('Percentage residuals for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
         percent_limits_2 = 20
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual [%]'
          plt.suptitle('Percentage residuals for remote-sensing reflectance: $R_{rs}$($\lambda$)')
          percent_limits = 16
          percent_limits_2 = 32
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual [%]'
          percent_limits = 40
          percent_limits_2 = 40
          plt.suptitle('Percentage residuals for normalized water-leaving radiance: $L_{wn}$($\lambda$)')

    xlab = 'Wavelength'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_PML, df_R ,bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, df_R ,bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    _resid_subplot(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, df_R ,bands)

   
    if spec_type == 'Lsky' :

        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
        
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits_2, df_NASA, df_R ,bands)
        
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits_2, df_CNR, df_R ,bands)
   
    elif spec_type == 'Lt' :
        
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
        
        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
    
                
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_CNR, df_R ,bands)
        
    elif spec_type == 'Rrs' :
          
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
        
        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
    
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_CNR, df_R ,bands)
   
    else: 

      subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
      index = 4
      _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
      
      subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
      index = 5
      _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_RBINS, df_R ,bands)
    
      subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
      index = 6
      _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_CNR, df_R ,bands)

        
    if spec_type == 'Ed' or  spec_type == 'nLw':
        subtitle  = 'NOAA: N = ' + str(np.sum(~np.isnan(df_NOAA['400'])))
        index = 7
        _resid_subplot(spec_type,subtitle, index, ylab, percent_limits_2, df_NOAA, df_R ,bands)
        
    
    plt.tight_layout()
    
    filename  =  path_output +  '/' + spec_type + '_resiudalsplot.png'
    plt.savefig(filename)
    
    return


def tabular_summary(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, df_CNR, df_NOAA, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to output tabular summary of results based on Tilstone 2020'''    
    
    # QC filtering
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
    df_CNR = df_CNR[Q_mask[Qtype]==1]
    df_NOAA = df_NOAA[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
        
    #  columns for table
    Institution = ['PML', 'NASA', 'TARTU', 'HEREON', 'RBINS', 'CNR', 'NOAA']
    
    Sensor = ['Seabird-HyperSAS', 'Seabird-HyperSAS', 'TriOS-RAMSES', 'TriOS-RAMSES', '', '', '']
    
    N_PML = np.sum(~np.isnan(df_PML['400']))
    N_NASA = np.sum(~np.isnan(df_NASA['400']))
    N_TARTU = np.sum(~np.isnan(df_TARTU['400']))
    N_HEREON = np.sum(~np.isnan(df_HEREON['400']))
    N_RBINS = np.sum(~np.isnan(df_RBINS['400']))
    N_CNR = np.sum(~np.isnan(df_CNR['400']))
    N_NOAA = np.sum(~np.isnan(df_NOAA['400']))
    N_meas = [N_PML, N_NASA, N_TARTU, N_HEREON, N_RBINS, N_CNR, N_NOAA]
    
    RMSD_442_PML = np.sqrt(np.nanmean((df_PML['442.5'] - df_R['442.5'])**2))
    RMSD_442_NASA = np.sqrt(np.nanmean((df_NASA['442.5'] - df_R['442.5'])**2))
    RMSD_442_TARTU =np.sqrt(np.nanmean((df_TARTU['442.5'] - df_R['442.5'])**2))
    RMSD_442_HEREON = np.sqrt(np.nanmean((df_HEREON['442.5'] - df_R['442.5'])**2))
    RMSD_442_RBINS = np.sqrt(np.nanmean((df_RBINS['442.5'] - df_R['442.5'])**2))
    RMSD_442 = [RMSD_442_PML, RMSD_442_NASA, RMSD_442_TARTU, RMSD_442_HEREON, RMSD_442_RBINS]
    
    RMSD_560_PML =  np.sqrt(np.nanmean((df_PML['560'] - df_R['560'])**2))
    RMSD_560_NASA =  np.sqrt(np.nanmean((df_NASA['560'] - df_R['560'])**2))
    RMSD_560_TARTU =  np.sqrt(np.nanmean((df_TARTU['560'] - df_R['560'])**2))
    RMSD_560_HEREON = np.sqrt(np.nanmean((df_HEREON['560'] - df_R['560'])**2))
    RMSD_560_RBINS = np.sqrt(np.nanmean((df_RBINS['560'] - df_R['560'])**2))
    RMSD_560 = [RMSD_560_PML, RMSD_560_NASA, RMSD_560_TARTU, RMSD_560_HEREON, RMSD_560_RBINS]
    
    RMSD_665_PML =  np.sqrt(np.nanmean((df_PML['665'] - df_R['665'])**2))
    RMSD_665_NASA =  np.sqrt(np.nanmean((df_NASA['665'] - df_R['665'])**2))
    RMSD_665_TARTU =  np.sqrt(np.nanmean((df_TARTU['665'] - df_R['665'])**2))
    RMSD_665_HEREON = np.sqrt(np.nanmean((df_HEREON['665'] - df_R['665'])**2))
    RMSD_665_RBINS = np.sqrt(np.nanmean((df_RBINS['665'] - df_R['665'])**2))
    RMSD_665 = [RMSD_665_PML, RMSD_665_NASA, RMSD_665_TARTU, RMSD_665_HEREON, RMSD_665_RBINS]
    
    RPD_442_PML =  100*np.nanmean((df_PML['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_NASA = 100*np.nanmean((df_NASA['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_TARTU =  100*np.nanmean((df_TARTU['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_HEREON = 100*np.nanmean((df_HEREON['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_RBINS = 100*np.nanmean((df_RBINS['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442 = [RPD_442_PML, RPD_442_NASA, RPD_442_TARTU, RPD_442_HEREON, RPD_442_RBINS]
    
    RPD_560_PML = 100*np.nanmean((df_PML['560'] - df_R['560'])/df_R['560'])
    RPD_560_NASA =  100*np.nanmean((df_NASA['560'] - df_R['560'])/df_R['560'])
    RPD_560_TARTU =  100*np.nanmean((df_TARTU['560'] - df_R['560'])/df_R['560'])
    RPD_560_HEREON = 100*np.nanmean((df_HEREON['560'] - df_R['560'])/df_R['560'])
    RPD_560_RBINS =  100*np.nanmean((df_RBINS['560'] - df_R['560'])/df_R['560'])
    RPD_560 = [RPD_560_PML, RPD_560_NASA, RPD_560_TARTU, RPD_560_HEREON, RPD_560_RBINS]
    
    RPD_665_PML =  100*np.nanmean((df_PML['665'] - df_R['665'])/df_R['665'])
    RPD_665_NASA =  100*np.nanmean((df_NASA['665'] - df_R['665'])/df_R['665'])
    RPD_665_TARTU =  100*np.nanmean((df_TARTU['665'] - df_R['665'])/df_R['665'])
    RPD_665_HEREON =  100*np.nanmean((df_HEREON['665'] - df_R['665'])/df_R['665'])
    RPD_665_RBINS =  100*np.nanmean((df_RBINS['665'] - df_R['665'])/df_R['665'])
    RPD_665 = [RPD_665_PML, RPD_665_NASA, RPD_665_TARTU, RPD_665_HEREON, RPD_665_RBINS]
    
    # convert to df format - 
    summary = pd.DataFrame() 
    summary['Institution'] = Institution 
    summary['Sensor type'] = Sensor
    summary['N'] = N_meas
    summary['RMSD 442.5'] = RMSD_442
    summary['RPD 442.5'] = RPD_442
    summary['RMSD 560'] = RMSD_560
    summary['RPD 560'] = RPD_560
    summary['RMSD 665'] = RMSD_665
    summary['RPD 665'] = RPD_665
      
    if spec_type != 'Rrs':
        summary['RMSD 442.5'] = summary['RMSD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 442.5'] = summary['RPD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RMSD 560'] = summary['RMSD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 560'] = summary['RPD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RMSD 665'] = summary['RMSD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 665'] = summary['RPD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))     
    elif spec_type == 'Rrs': 
       summary['RPD 442.5'] = summary['RPD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
       summary['RPD 560'] = summary['RPD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
       summary['RPD 665'] = summary['RPD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))

    filename  =  path_output +  '/' + spec_type + '_summary.csv'
    summary.to_csv(filename, na_rep ='NaN', index = False)
    filename2  =  path_output +  '/' + spec_type + '_summary.png'

    dfi.export(summary.style.hide(axis='index'), filename2)
    
    return summary

def filter_by_azimuth(df_PML, df_NASA, df_RBINS, tol=1):
    '''function to filter NASA and RBINS via azimuth'''   
    
    delta_phi_NASA = abs(df_NASA['azimuth']) - abs(df_PML['azimuth'])
    delta_phi_RBINS = abs(df_RBINS['azimuth']) - abs(df_PML['azimuth'])
   
    # delta_phi_NASA = df_NASA['azimuth'] - df_PML['azimuth']
    # delta_phi_RBINS = df_RBINS['azimuth'] - df_PML['azimuth']
    
    for i in range(78):
       if abs(delta_phi_NASA[i]) > tol:
         df_NASA.iloc[i] = np.nan
       if abs(delta_phi_RBINS[i]) > tol:
         df_RBINS.iloc[i] = np.nan
         
    return df_NASA, df_RBINS

def azimuth_plot():
    
    plt.figure()
    plt.title('Relative azimuth (absolute value)')
    plt.plot(Ed_PML.index, abs(Ed_PML['azimuth']), label= 'PML, TARTU, HEREON',color='orange',linewidth = 3)
    plt.plot(Ed_NASA.index, abs(Ed_NASA['azimuth']), label= 'NASA', linestyle='dashed',color='blue',linewidth = 3)
    plt.plot(Ed_RBINS.index, abs(Ed_RBINS['azimuth']), label= 'RBINS', linestyle='dashed', color='red',linewidth =3)
    plt.xlabel('Station number')
    plt.ylabel('|$\phi$| (deg)')
    plt.legend()
    
    #plt.figure()
    #plt.title('Relative azimuth')
    #plt.plot(Ed_PML.index, Ed_PML['azimuth'], label= 'PML, TARTU, HEREON',color='orange',  linewidth = 3)
    #plt.plot(Ed_NASA.index, Ed_NASA['azimuth'], label= 'NASA', linestyle='dashed',color='blue', linewidth = 3)
    #plt.plot(Ed_RBINS.index, Ed_RBINS['azimuth'], label= 'RBINS', linestyle='dashed', color='red', linewidth = 3)
    #plt.xlabel('Station number')
    #plt.ylabel('$\phi$ (deg)')
    #plt.legend()
    
    filename  =  path_output +  '/azimuth.png'
    plt.savefig(filename)

    return


if __name__ == '__main__':
    
    # options
    # baseline =
    # QC method =
    
    # path to data + output
    dir_data = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/DataSubmissions'
    path_output = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Output'

    # team submissions
    path_PML = dir_data + '/PML/FICE_submission_V2/FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_L1timesWindspeeds_phi_output.csv'
    path_NASA = dir_data + '/NASA/FRM4SOC_2_FICE_22_AAOT_pySAS_Sentinel3A_rev1.csv'
    path_TARTU = dir_data + '/UniTartu/averages.csv'  
    path_RBINS = dir_data + '/RBINS/FRM4SOC2_FICE_RBINS_2022-09-16.csv'
    path_HEREON = dir_data + '/HEREON/FRM4SOC2_FICE22_HEREON_DATA_OLCI'
    path_NOAA = dir_data + '/NOAA/NOAA_Hyperpro_sheet2.csv'
    path_CNR = dir_data +  '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_averages+stdev.csv'

    # addtional data (QC + references)    
    path_QC = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Aeronet_QC_mask/FRM4SOC-AAOT_V3_ope.txt'
    path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference'    

    #OLCI bands
    bands = [str(400), str(412.5), str(442.5),	str(490), str(510), str(560), str(620),	str(665), str(673.75), str(681.25), str(708.75), str(753.75), str(761.25), str(764.375), str(767.5), str(778.75), str(865), str(885), str(900)]
 
    # Read data - each team has own file reader to homogenize data
    Ed_PML, Lsky_PML, Lt_PML, Rrs_PML, Rrs_std_PML, nLw_PML = read_PML_data(path_PML, bands)
    Ed_NASA, Lsky_NASA, Lt_NASA, Rrs_NASA, Rrs_std_NASA, nLw_NASA = read_NASA_data(path_NASA, bands)
    Ed_TARTU, Lsky_TARTU, Lt_TARTU, Rrs_TARTU, Rrs_std_TARTU, nLw_TARTU = read_TARTU_data(path_TARTU, bands)
    Ed_HEREON, Lsky_HEREON, Lt_HEREON, Rrs_HEREON, Rrs_std_HEREON, nLw_HEREON = read_HEREON_data(path_HEREON, bands)
    Ed_RBINS, Lsky_RBINS, Lt_RBINS, Rrs_RBINS, Rrs_std_RBINS, nLw_RBINS = read_RBINS_data(path_RBINS, bands)
    Ed_NOAA, Lsky_NOAA, Lt_NOAA, Rrs_NOAA, nLw_NOAA  = read_NOAA_data(path_NOAA, bands, Ed_PML) # PML timestamps used to reference staion no.
    Ed_CNR, Lsky_CNR, Lt_CNR, Rrs_CNR, Rrs_std_CNR, nLw_CNR = read_CNR_data(path_CNR, bands)
    
    # Azimuth filtering for NASA and RBINS
    Lsky_NASA, Lsky_RBINS = filter_by_azimuth(Lsky_PML, Lsky_NASA, Lsky_RBINS)
    Lt_NASA, Lt_RBINS = filter_by_azimuth(Lt_PML, Lt_NASA, Lt_RBINS)
    Rrs_NASA, Rrs_RBINS = filter_by_azimuth(Rrs_PML, Rrs_NASA, Rrs_RBINS)
    
    # References/ baselines
    # Ed_R = baseline_average('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON)
    # Lsky_R = baseline_average('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON)
    # Lt_R = baseline_average('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON)
    #Rrs_R = baseline_average('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON)
    # Rrs_std_R = baseline_average('Rrs', Rrs_std_PML, Rrs_std_NASA, Rrs_std_TARTU, Rrs_std_HEREON)  
   
    Ed_R = baseline_average_V2('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON)
    Lsky_R = baseline_average_V2('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON)
    Lt_R = baseline_average_V2('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON)
    Rrs_R = baseline_average_V2('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON)
    Rrs_std_R = baseline_average_V2('Rrs', Rrs_std_PML, Rrs_std_NASA, Rrs_std_TARTU, Rrs_std_HEREON)  
    nLw_R = read_Aeronet_nLw(path_NLW, Ed_R, bands)    
    
    # Quality control
    Q_mask = QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML) 
    plot_dataandQCmasks()
    azimuth_plot()
    
    # results plots    
    plot_scatter('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Ed_CNR, Ed_NOAA, Q_mask, Qtype = 'QC_CI_3')   
    plot_scatter('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_scatter('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_scatter('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS,  Rrs_CNR, Rrs_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_scatter('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, Q_mask, Qtype = 'QC_CI_3') 
        
    plot_residuals('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Ed_CNR, Ed_NOAA, Q_mask, Qtype = 'QC_CI_3')    
    plot_residuals('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_residuals('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_residuals('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, Rrs_CNR, Rrs_NOAA, Q_mask, Qtype = 'QC_CI_3') 
    plot_residuals('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, Q_mask, Qtype = 'QC_CI_3') 
 
    # Ed_table = tabular_summary('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Q_mask, Qtype = 'QC_AOC_3')    
    # Lsky_table = tabular_summary('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Q_mask, Qtype = 'QC_AOC_3') 
    # Lt_table = tabular_summary('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Q_mask, Qtype = 'QC_AOC_3') 
    # Rrs_table = tabular_summary('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, Q_mask, Qtype = 'QC_AOC_3') 
    # nLw = tabular_summary('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, Q_mask, Qtype = 'QC_AOC_3') 
    
    # Time series for 19/0