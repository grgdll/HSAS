#!/usr/bin/env python

'''

Code to extract matchup Rrs from AMT hypersas system

Author: hek
history: hek created for bob 22/03/18
	 hek edited August 2018 for consistency with AMT26 and for AMT4SFRM matchups 
version: 2
'''

import csv
import numpy as np
import os, sys, shutil
import fnmatch
import datetime
import scipy.io as io
import matplotlib.pyplot as plt
from matplotlib import dates

### Functions

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return [idx]

### Set directories


in_dir_H='/data/datasets/cruise_data/active/AMT27/Public_Read_Only_Copy/AMT27_DY084/Optics_group/Processed/rad_processing/data/Hsas/Processed_final/L2/'
out_file_H='/data/datasets/cruise_data/active/AMT27/Public_Read_Only_Copy/AMT27_DY084/Optics_group/Processed/rad_processing/data/Hsas/Processed_final/matchup_extractions/AMT27_HSAS_all_v3.txt'
Lambda=np.array(range(400,702,2))


# remove file if it exists
if os.path.exists(out_file_H):
    os.remove(out_file_H)

# open new file
f = open(out_file_H,'w')

# create header
HeaderH = 'latitude,longitude,Hsas_time'
for l in np.arange(0,len(Lambda)):
    HeaderH = HeaderH+', RRS ['+str(Lambda[l])+' nm]'

# write header
f.write(HeaderH+'\n')

years = [2017]
days  = np.arange(267,309)

for year in years:
   for day in days:
       if day == 0:
           continue
       print(year,' : ',day)
       mydate = datetime.datetime.strptime(str(year)+str(day),'%Y%j')
       YYYY = str(mydate.year)
       MM = str(mydate.month).zfill(2)
       DD = str(mydate.day).zfill(2)
       JJ = str(day).zfill(3)
       
       data_file = in_dir_H+YYYY+'_'+MM+'_'+DD+'_Hsas_L2_v1_Cont.mat'
       print data_file
       
       if os.path.exists(data_file):
          print('Opening: '+data_file)
       else:
           continue
       print('Opening: '+data_file)
       
       load_file = io.loadmat(data_file,squeeze_me=True, struct_as_record=False)
       Hsas_time = load_file['L2'].time
       Hsas_data = load_file['L2'].Rrs.data[:,25:176]
       Hsas_lat = load_file['L2'].gps.lat
       Hsas_lon = load_file['L2'].gps.lon

       #Hsas_dt=np.zeros(np.shape(Hsas_time))
       Hsas_dt=[]
           
       for i in range(0,len(Hsas_time)):
          Hsas_datetime=(datetime.datetime.fromordinal(int(Hsas_time[i]))+datetime.timedelta(seconds=86400*(Hsas_time[i]-float(int(Hsas_time[i]))))-datetime.timedelta(days=366))
          Hsas_dt.append(Hsas_datetime.strftime("%Y-%m-%dT%H:%M:%SZ"))
          # concat vars and write
          lon  = str(Hsas_lon[i])
          lat  = str(Hsas_lat[i])
          dt   = str(Hsas_dt[i])
          data = [str(x) for x in Hsas_data[i,:]]
          row = lon+','+lat+','+dt+','+','.join(data)
          f.write(row+'\n')

f.close()

