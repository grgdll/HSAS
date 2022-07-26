#!/usr/bin/env python

'''

Code to calculate uncertainties for underway Rrs and write out to single text file

Author: hek
history: hek editted to calculate uncertainty for AMT27 (different calibration process from 26)
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
import pandas

### Functions
########## Define functions for individual uncertainty calculations ############

def env_unc(meas_counts):
    env_term_out=np.zeros((np.shape(meas_counts)))
    for i in np.arange(0,np.shape(meas_counts)[1]):
        meas_slice=meas_counts[:,i]
        env_term=pandas.rolling_std(meas_slice,3)**2
        env_term2 = pandas.rolling_std(meas_slice[::-1],3)**2
        env_term2 = env_term2[::-1]
        env_term[:2]=env_term2[:2]
        env_term_out[:,i]=env_term
    return env_term_out
	
def ins_meas_unc(sig_calcoeff_sq, meas_counts, sig_env, ins):
    count_term=(sig_calcoeff_sq)
    coeff_term=(sig_calcoeff_sq)
    sig_ins_sq = count_term+coeff_term+sig_env

    if ins == 'es':
        sl_term=(meas_counts*0.025)**2
        cos_term=(meas_counts*0.02)**2
        pol_term=(meas_counts*0.01)**2
        sig_ins_sq=sig_ins_sq+sl_term+cos_term+pol_term
        pc_meas_var=count_term/sig_ins_sq
        pc_cal_coef=coeff_term/sig_ins_sq
        pc_ct=cos_term/sig_ins_sq
        pc_sl=sl_term/sig_ins_sq
        pc_pt=pol_term/sig_ins_sq
        pc_env=sig_env/sig_ins_sq
    else:
        sl_term=(meas_counts*0.025)**2
        cos_term=0
        pol_term=(meas_counts*0.01)**2
        sig_ins_sq=sig_ins_sq+sl_term+cos_term+pol_term
        pc_meas_var=count_term/sig_ins_sq
        pc_cal_coef=coeff_term/sig_ins_sq
        pc_ct=cos_term/sig_ins_sq
        pc_sl=sl_term/sig_ins_sq
        pc_pt=pol_term/sig_ins_sq
        pc_env=sig_env/sig_ins_sq

    
    return sig_ins_sq, pc_meas_var, pc_cal_coef, pc_ct, pc_sl, pc_pt, pc_env
                    
                    
def calc_cc_terms(lt_in, li_in, es_in):
    cc_lt_li=np.zeros((np.shape(lt_in)[1]))
    cc_lt_es=np.zeros((np.shape(lt_in)[1]))
    cc_li_es=np.zeros((np.shape(lt_in)[1]))
    for i in np.arange(0,np.shape(lt_in)[1]):
        cc_lt_li[i] = np.corrcoef(lt_in[:,i],li_in[:,i])[0,1]
        cc_lt_es[i] = np.corrcoef(lt_in[:,i],es_in[:,i])[0,1]
        cc_li_es[i] = np.corrcoef(li_in[:,i],es_in[:,i])[0,1]
    return cc_lt_li, cc_lt_es, cc_li_es

def calc_combined(es_samp, sig_es, lt_samp, sig_lt, li_samp, sig_li, rho_samp, cc_lt_li, cc_lt_es, cc_li_es, rrs_in):
    rrs_lt_term = ((1/es_samp)*sig_lt)**2
    rrs_li_term = ((-rho_samp/es_samp)*sig_li)**2
    rrs_es_term = (((-lt_samp+rho_samp*li_samp)/es_samp**2)*sig_es)**2
    rrs_rho_term = ((-li_samp/es_samp)*(rho_samp*0.025))**2
    c_lt_li_term = 2*(1/es_samp)*(-rho_samp/es_samp)*sig_lt*sig_li*cc_lt_li
    c_lt_es_term = 2*(1/es_samp)*((-lt_samp+rho_samp*li_samp)/es_samp**2)*sig_lt*sig_es*cc_lt_es
    c_li_es_term = 2*((-rho_samp/es_samp)*((-lt_samp+rho_samp*li_samp)/es_samp**2))*sig_li*sig_es*cc_li_es
    
    #sig_rrs_nocorr=rrs_lt_term + rrs_li_term + rrs_es_term + rrs_rho_term
    
    sig_rrs_sq=rrs_lt_term + rrs_li_term + rrs_es_term + rrs_rho_term+c_lt_li_term + c_lt_es_term + c_li_es_term
    sig_rrs_pie=rrs_lt_term + rrs_li_term + rrs_es_term + rrs_rho_term+abs(c_lt_li_term) + abs(c_lt_es_term) + abs(c_li_es_term)

    pc_lt=rrs_lt_term/sig_rrs_pie
    pc_li=rrs_li_term/sig_rrs_pie
    pc_es=rrs_es_term/sig_rrs_pie
    pc_rho=rrs_rho_term/sig_rrs_pie
    pc_cltli=abs(c_lt_li_term)/sig_rrs_pie
    pc_cltes=abs(c_lt_es_term)/sig_rrs_pie
    pc_clies=abs(c_li_es_term)/sig_rrs_pie

    return sig_rrs_sq, pc_lt, pc_li, pc_es, pc_rho,pc_cltli,pc_cltes,pc_clies


def perform_unc_summation(Hsas_wv_new, t_in):
    ##### Run through each step of uncertainty derivation #####

   # 1. Calculate the uncertainty in the calibration coefficient 

    t_lt_unc_sq=(t_in['L2'].instr.Lt.data[:,25:176]*0.01)**2
    t_li_unc_sq=(t_in['L2'].instr.Li.data[:,25:176]*0.01)**2
    t_es_unc_sq=(t_in['L2'].instr.Es.data[:,25:176]*0.01)**2

    print np.shape

   # 2. Calculate the measurement uncertainty for each measurement (lt, li, es)

    # Load dummy variable to get size correct for initialising arrays

    shape_init=t_in['L2'].instr.Es.data[:,25:176]

    t_sig_lt=np.zeros((np.shape(shape_init)))
    t_sig_li=np.zeros((np.shape(shape_init)))
    t_sig_es=np.zeros((np.shape(shape_init)))

    pc_meas_lt=np.zeros((np.shape(shape_init)))
    pc_dark_lt=np.zeros((np.shape(shape_init)))
    pc_cal_lt=np.zeros((np.shape(shape_init)))
    pc_ct_lt=np.zeros((np.shape(shape_init)))
    pc_sl_lt=np.zeros((np.shape(shape_init)))
    pc_pt_lt=np.zeros((np.shape(shape_init)))
    pc_env_lt=np.zeros((np.shape(shape_init)))

    pc_meas_li=np.zeros((np.shape(shape_init)))
    pc_dark_li=np.zeros((np.shape(shape_init)))
    pc_cal_li=np.zeros((np.shape(shape_init)))
    pc_ct_li=np.zeros((np.shape(shape_init)))
    pc_sl_li=np.zeros((np.shape(shape_init)))
    pc_pt_li=np.zeros((np.shape(shape_init)))
    pc_env_li=np.zeros((np.shape(shape_init)))

    pc_meas_es=np.zeros((np.shape(shape_init)))
    pc_dark_es=np.zeros((np.shape(shape_init)))
    pc_cal_es=np.zeros((np.shape(shape_init)))
    pc_ct_es=np.zeros((np.shape(shape_init)))
    pc_sl_es=np.zeros((np.shape(shape_init)))
    pc_pt_es=np.zeros((np.shape(shape_init)))
    pc_env_es=np.zeros((np.shape(shape_init)))


    t_es=t_in['L2'].instr.Es.data[:,25:176]
    t_lt=t_in['L2'].instr.Lt.data[:,25:176]
    t_li=t_in['L2'].instr.Li.data[:,25:176]
    t_time=t_in['L2'].time

    # Calculate uncertainty associated with environmental factors (based on moving std with time)

    env_term_lt=env_unc(t_lt)
    env_term_li=env_unc(t_li)
    env_term_es=env_unc(t_es)

    for i in range(0,len(t_time)):
        t_sig_lt[i,:], pc_meas_lt[i,:], pc_cal_lt[i,:], pc_ct_lt[i,:], pc_sl_lt[i,:], pc_pt_lt[i,:], pc_env_lt[i,:] = ins_meas_unc(t_lt_unc_sq[i,:], t_in['L2'].instr.Lt.data[i,25:176], env_term_lt[i,:],'lt')
        t_sig_li[i,:], pc_meas_li[i,:], pc_cal_li[i,:], pc_ct_li[i,:], pc_sl_li[i,:], pc_pt_li[i,:], pc_env_li[i,:] = ins_meas_unc(t_li_unc_sq[i,:], t_in['L2'].instr.Li.data[i,25:176], env_term_li[i,:],'li')
        t_sig_es[i,:], pc_meas_es[i,:], pc_cal_es[i,:], pc_ct_es[i,:], pc_sl_es[i,:], pc_pt_es[i,:], pc_env_es[i,:] = ins_meas_unc(t_es_unc_sq[i,:], t_in['L2'].instr.Es.data[i,25:176], env_term_es[i,:],'es')

    t_sig_lt=np.sqrt(t_sig_lt)
    t_sig_li=np.sqrt(t_sig_li)
    t_sig_es=np.sqrt(t_sig_es)


    # 3. Calculate the uncertainty in the resulting Rrs

    # Determine correlation coefficients for covarying terms (i.e. all light based measurements!!)

    cc_lt_li, cc_lt_es, cc_li_es=calc_cc_terms(t_in['L2'].instr.Lt.data[:,25:176], t_in['L2'].instr.Li.data[:,25:176], t_in['L2'].instr.Es.data[:,25:176])

    # Perform summation of uncertainties

    t_sig_rrs_sq=np.zeros((np.shape(shape_init)))
    t_sig_rrs=np.zeros((np.shape(shape_init)))
    pc_lt=np.zeros((np.shape(shape_init)))
    pc_li=np.zeros((np.shape(shape_init)))
    pc_es=np.zeros((np.shape(shape_init)))
    pc_rho=np.zeros((np.shape(shape_init)))
    pc_cltli=np.zeros((np.shape(shape_init)))
    pc_cltes=np.zeros((np.shape(shape_init)))
    pc_clies=np.zeros((np.shape(shape_init)))

    for i in range(0,len(t_time)):
        t_sig_rrs_sq[i,:],pc_lt[i,:], pc_li[i,:], pc_es[i,:], pc_rho[i,:],pc_cltli[i,:],pc_cltes[i,:],pc_clies[i,:]=calc_combined(t_in['L2'].instr.Es.data[i,25:176], t_sig_es[i,:], t_in['L2'].instr.Lt.data[i,25:176], t_sig_lt[i,:], t_in['L2'].instr.Li.data[i,25:176], t_sig_li[i,:], t_in['L2'].rho_fitted[i,0], cc_lt_li[:], cc_lt_es[:], cc_li_es[:], t_in['L2'].Rrs.data[i,25:176])

    # Calculate uncertainties as percent - OR IS THIS CORRECT? I THINK THIS IS (error comes out in units, so need to calculate as a percentage of each measurement?

    unc=(np.sqrt(t_sig_rrs_sq)/t_in['L2'].Rrs.data[:,25:176])*100

    #mean_perwv=np.mean(unc[:,0:101],0)
    mean_persamp=np.mean(unc,1)

    return  mean_persamp, unc


# Main program

### Set directories

in_dir_H='/data/datasets/cruise_data/active/AMT24/DallOlmo/HSAS/Processed_final/L2/'

Lambda=np.array(range(400,702,2))

out_file_H='/data/datasets/cruise_data/active/AMT24/DallOlmo/HSAS/Processed_final/uncertainty_estimates/amt24_uncertainty_estimates.txt'

# remove output file if it exists
if os.path.exists(out_file_H):
    os.remove(out_file_H)

# open new file
f = open(out_file_H,'w')

# create header
Headerbase = 'latitude,longitude,time,totunc'
for l in np.arange(0,len(Lambda)):
    Headerbase = Headerbase+', RRSUNC['+str(Lambda[l])+' nm]'

# write header
f.write(Headerbase+'\n')

years = [2014]
days  = np.arange(268,305)

Hsas_tot_unc_all=[]

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
       
       if os.path.exists(data_file):
          print('Opening: '+data_file)
       else:
           continue
       
       Hsas_file = io.loadmat(data_file,squeeze_me=True, struct_as_record=False)
       Hsas_time = Hsas_file['L2'].time
       Hsas_data = Hsas_file['L2'].Rrs.data[:,25:176]
       Hsas_lat = Hsas_file['L2'].gps.lat
       Hsas_lon = Hsas_file['L2'].gps.lon
       
       # Deal with wavelength for Hsas

       Hsas_wv = Hsas_file['L2'].wv
       Hsas_wv_new = Hsas_wv[25:176]

       #Hsas_dt=np.zeros(np.shape(Hsas_time))
       Hsas_dt=[]
       
       #Pass data to calibration function
       
       Hsas_tot_unc, Hsas_rrs_unc = perform_unc_summation(Hsas_wv_new, Hsas_file)
       
       Hsas_tot_unc_all.append(Hsas_tot_unc)
       for i in range(0,len(Hsas_time)):
          Hsas_datetime=(datetime.datetime.fromordinal(int(Hsas_time[i]))+datetime.timedelta(seconds=86400*(Hsas_time[i]-float(int(Hsas_time[i]))))-datetime.timedelta(days=366))
          Hsas_dt.append(Hsas_datetime.strftime("%Y-%m-%dT%H:%M:%SZ"))
          # concat vars and write
          lon  = str(Hsas_lon[i])
          lat  = str(Hsas_lat[i])
          dt   = str(Hsas_dt[i])
          tot  = str(Hsas_tot_unc[i])
          data = [str(x) for x in Hsas_rrs_unc[i,:]]
          row = lon+','+lat+','+dt+','+tot+','+','.join(data)
          f.write(row+'\n')

f.close()

