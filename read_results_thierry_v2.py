#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 04:30:22 2017

@author: jmilli
"""


import os
import numpy as np
import irdisDataHandler as i
import pdb
import vip
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
from scipy.interpolate import interp1d

local = False

if local:
    pathRoot = '/diskb/Data'
else:
    pathRoot = '/Volumes/SHARDDS_data'

#%%
# We read the measured strehl from Thierry
path = os.path.join(pathRoot,'PSF_analysis_final/result_analysis_thierry')

#name = fits.getdata(os.path.join(path,'name.fits'))
#hdulist_name = fits.open(os.path.join(path,'name.fits'))
#hdulist_name.info()
target_names = ['AG_Tri_2nd_epoch','AG_Tri_3rd_epoch','AG_Tri_4th_epoch', 'AG_Tri', \
        'FomC_2nd_epoch','FomC_3rd_epoch','bad_FomC_4th_epoch','FomC', 'HD105', \
        'HD203', 'HD377', 'HD3670', 'HD9672', 'HD10472_2nd_epoch', 'HD10472', \
        'HD10638', 'HD13246', 'HD14082B_2nd_epoch','HD14082B', 'HD15257', \
        'HD16743', 'HD17390', 'HD22179_2nd_epoch', 'HD22179', 'HD24636', 'HD25457', \
        'HD31392', 'HD35650', 'HD37484_2nd_epoch', 'HD37484', 'HD38206_2nd_epoch', \
        'HD38206', 'HD38206' , 'HD40540' , 'HD53842' , 'HD60491' , 'HD69830' , \
        'HD71722' , 'HD73350' , 'HD76582' , 'HD80950' , 'HD82943_2nd_epoch' ,\
        'HD82943_3rd_epoch', 'HD82943_4th_epoch', 'HD82943','HD84075' , 'HD107649' , \
        'HD114082' , 'HD120534_2nd_epoch' , 'HD120534', 'HD122652_2nd_epoch' , \
        'HD133803_2nd_epoch' , 'HD133803', 'HD135599_3rd_epoch' , 'HD138965' , \
        'HD145229' , 'HD157728' , 'HD164249A' , 'HD172555' , 'HD181296' , \
        'HD182681' , 'HD192758_2nd_epoch' , 'HD192758', 'HD201219' , \
        'HD205674_2nd_epoch' , 'HD205674', 'HD206893', 'HD218340', 'HD221853', \
        'HD274255', 'HIP63942']

print(target_names)

nb_targets = len(target_names)
print(nb_targets)

psfall_norm_at_SR = fits.getdata(os.path.join(path,'psfall_norm_at_SR.fits'))

sr = fits.getdata(os.path.join(path,'sr.fits'))
print(len(sr)) #71

ds9=vip.fits.vipDS9()
#ds9.display(psfall_norm_at_SR)
ds9.display(psfall_norm_at_SR[sr>0.8],psfall_norm_at_SR[sr<0.6],psfall_norm_at_SR[sr<0.6])

print('\n'.join(np.asarray(target_names)[sr>0.8][5:]))
print(np.asarray(sr)[sr>0.8][5:])

print(np.asarray(target_names)[sr<0.6])
print(np.asarray(sr)[sr<0.6])

# illustration
# 2 for LWE


#%%
# We read the estimation from sparta along with other keywords, both for the PSF 
# and the O frames.

# Creation of one overall dictionary to store the mean, median, rms, max and min values
complete_dico = {'target_name':target_names,'measured_strehl':sr}

# We instantiate first all available fields:
for target_name in target_names:
    path_target = os.path.join(pathRoot,os.path.join('survey_disk',target_name))
    path_pipeline = os.path.join(path_target,'pipeline')
    file_sparta_stats_F = os.path.join(path_pipeline,target_name+'_keywords_statistics_F.csv')
    pd_target_sparta_stats_F = pd.read_csv(file_sparta_stats_F)
    # We add all the fields from the keywords_statistics_F.csv file
    fields = list(pd_target_sparta_stats_F['name'])
    for field in fields:
        if field+'_F' not in complete_dico.keys():
            complete_dico[field+'_F'] = np.ones(nb_targets)*np.nan
    # We add all the fields from the keywords_statistics_O.csv file
    file_sparta_stats_O = os.path.join(path_pipeline,target_name+'_keywords_statistics_O.csv')
    if os.path.isfile(file_sparta_stats_O):
        pd_target_sparta_stats_O = pd.read_csv(file_sparta_stats_O)    
        fields = list(pd_target_sparta_stats_O['name'])
        for field in fields:
            fields_to_check = [field+'_O_mean',field+'_O_median',field+'_O_rms',\
                           field+'_O_max',field+'_O_min']
            for field_to_check in fields_to_check:
                if field_to_check not in complete_dico.keys():
                    complete_dico[field_to_check] = np.ones(nb_targets)*np.nan
    # We add all the fields from the simbad.csv file
    file_simbad = os.path.join(path_pipeline,target_name+'_simbad_info.csv')
    if os.path.isfile(file_simbad):
        pd_simbad = pd.read_csv(file_simbad)
        for field in pd_simbad.keys():            
            if field not in complete_dico.keys() and field.startswith('simbad_FLUX'):
                complete_dico[field] = np.ones(nb_targets)*np.nan
   
# We populate this dictionary
for itarget,target_name in enumerate(target_names):
    print(target_name)
    path_target = os.path.join(pathRoot,os.path.join('survey_disk',target_name))
    path_pipeline = os.path.join(path_target,'pipeline')
    file_sparta_stats_F = target_name+'_keywords_statistics_F.csv'
    pd_target_sparta_stats_F = pd.read_csv(os.path.join(path_pipeline,\
                                                      file_sparta_stats_F))
    for ikeyword,keyword in enumerate(pd_target_sparta_stats_F['name']):
        complete_dico[keyword+'_F'][itarget] = pd_target_sparta_stats_F['mean'][ikeyword]
    file_sparta_stats_O = os.path.join(path_pipeline,target_name+'_keywords_statistics_O.csv')
    if os.path.isfile(file_sparta_stats_O):
        pd_target_sparta_stats_O = pd.read_csv(file_sparta_stats_O)
        for ikeyword,keyword in enumerate(pd_target_sparta_stats_O['name']):
            complete_dico[keyword+'_O_mean'][itarget] = pd_target_sparta_stats_O['mean'][ikeyword]
            complete_dico[keyword+'_O_median'][itarget] = pd_target_sparta_stats_O['median'][ikeyword]
            complete_dico[keyword+'_O_rms'][itarget] = pd_target_sparta_stats_O['rms'][ikeyword]
            complete_dico[keyword+'_O_max'][itarget] = pd_target_sparta_stats_O['max'][ikeyword]
            complete_dico[keyword+'_O_min'][itarget] = pd_target_sparta_stats_O['min'][ikeyword]
    else:
        print('{0:s} does not exist'.format(file_sparta_stats_O))
    file_simbad = os.path.join(path_pipeline,target_name+'_simbad_info.csv')
    if os.path.isfile(file_simbad):
        pd_simbad = pd.read_csv(file_simbad)
        for field in pd_simbad.keys():
            if field.startswith('simbad_FLUX'):
                complete_dico[field][itarget] = str(pd_simbad[field][0])
    else:
        print('{0:s} does not exist'.format(file_simbad))

# Conversion to pandas
pd_total = pd.DataFrame(complete_dico)

# We save everything
pd_total.to_csv(os.path.join(path,'complete_data.csv'),index=False)


#%%
# We open the final csv summary files

pd_total = pd.read_csv(os.path.join(path,'complete_data.csv'))

#%%
# Analysis

#s=np.power(10,-(pd_total['simbad_FLUX_V']-12.5)/2.5)
#s=1000/pd_total['simbad_FLUX_V']

plt.scatter(pd_total['measured_strehl'],complete_dico['strehl_F'],\
            c=pd_total['simbad_FLUX_V'],s=80,cmap=plt.cm.coolwarm)# cmap='gray')
plt.plot([0.2,1],[0.2,1],'--',color='black')
plt.xlabel('Measured Strehl')
plt.ylabel('RTC Strehl')
plt.axis([0.2,1,0.2,1])
#plt.axis([0,1,0,1])
#plt.axis('equal')
#plt.grid()
#plt.colorbar()
cbar = plt.colorbar()
cbar.set_label('V magnitude', rotation=270)
plt.savefig(os.path.join(path,'comparison_sparta_RTC.pdf'))


#%% We read the contrast files.


keys_to_add = ['contrast_cadi_200mas','contrast_cadi_500mas','contrast_cadi_1000mas',\
               'contrast_pca_200mas','contrast_pca_500mas','contrast_pca_1000mas',\
               'contrast_raw_coro_200mas','contrast_raw_coro_500mas','contrast_raw_coro_1000mas']
for k in keys_to_add:
    pd_total[k] = np.ones(nb_targets)*np.nan

px=0.01225

for itarget,target_name in enumerate(target_names):
    print(target_name)
    path_target = os.path.join(pathRoot,os.path.join('survey_disk',target_name))
    path_pipeline = os.path.join(path_target,'pipeline')
    file_contrast_cadi = os.path.join(path_pipeline,\
                        target_name+'_contrast_curve_199x199_left_sorted_cadi.csv')
    if os.path.isfile(file_contrast_cadi):
        pd_target_contrast_cadi = pd.read_csv(file_contrast_cadi)
        contrast = pd_target_contrast_cadi['sensitivity (Student)']
        sep_arcsec = pd_target_contrast_cadi['distance']*px
        contrast_function = interp1d(sep_arcsec,contrast,kind='linear',bounds_error=False,fill_value=np.nan)
        pd_total['contrast_cadi_200mas'][itarget] = contrast_function(0.2) 
        pd_total['contrast_cadi_500mas'][itarget] = contrast_function(0.5) 
        pd_total['contrast_cadi_1000mas'][itarget] = contrast_function(1.0) 
    else:
        print('No contrast file found: {0:s}'.format(file_contrast_cadi))
    file_contrast_pca = os.path.join(path_pipeline,\
                        target_name+'_contrast_curve_199x199_left_sorted_pca_klip_010.csv')
    if os.path.isfile(file_contrast_pca):
        pd_target_contrast_pca = pd.read_csv(file_contrast_pca)
        contrast = pd_target_contrast_pca['sensitivity (Student)']
        sep_arcsec = pd_target_contrast_pca['distance']*px
        contrast_function = interp1d(sep_arcsec,contrast,kind='linear',bounds_error=False,fill_value=np.nan)
        pd_total['contrast_pca_200mas'][itarget] = contrast_function(0.2) 
        pd_total['contrast_pca_500mas'][itarget] = contrast_function(0.5) 
        pd_total['contrast_pca_1000mas'][itarget] = contrast_function(1.0) 
    else:
        print('No contrast file found: {0:s}'.format(file_contrast_pca))
    file_contrast_raw = os.path.join(path_pipeline,\
                        target_name+'_contrast_curve_199x199_left_sorted_raw_coro.csv')
    if os.path.isfile(file_contrast_raw):
        pd_target_contrast_pca = pd.read_csv(file_contrast_raw)
        contrast = pd_target_contrast_pca['sensitivity (Student)']
        sep_arcsec = pd_target_contrast_pca['distance']*px
        contrast_function = interp1d(sep_arcsec,contrast,kind='linear',bounds_error=False,fill_value=np.nan)
        pd_total['contrast_raw_coro_200mas'][itarget] = contrast_function(0.2) 
        pd_total['contrast_raw_coro_500mas'][itarget] = contrast_function(0.5) 
        pd_total['contrast_raw_coro_1000mas'][itarget] = contrast_function(1.0) 
    else:
        print('No contrast file found: {0:s}'.format(file_contrast_raw))
        

#pd_contrast = pd.DataFrame(dico_contrast)
pd_total.to_csv(os.path.join(path,'complete_complete_data.csv'),index=False)

pd_total = pd.read_csv(os.path.join(path,'complete_complete_data.csv'))

plt.plot(pd_total['simbad_FLUX_V'],pd_total['strehl_F'],'o',color='rosybrown',label='RTC')
plt.plot(pd_total['simbad_FLUX_V'],pd_total['measured_strehl'],'o',color='black',label='measured')
plt.xlabel('V magnitude')
plt.ylabel('Strehl')
plt.grid()
plt.legend()
plt.axis([4,13,0.3,1])
plt.savefig(os.path.join(path,'SHARDDS_comparison_strehl_Vmag.pdf'))


plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_1000mas'],'go',label='1000mas')
plt.xlabel('Measured Strehl')
plt.ylabel('$5\sigma$ post-ADI (classical) contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_contrast_vs_strehl_cadi.pdf'))

# we add a trend line
sr= pd_total['measured_strehl']
c_200 = pd_total['contrast_cadi_200mas']
c_500 = pd_total['contrast_cadi_500mas']
c_1000 = pd_total['contrast_cadi_1000mas']
good =np.isfinite(c_200) 
c_200 = np.log10(c_200[good])
c_500 = np.log10(c_500[good])
c_1000 = np.log10(c_1000[good])
sr=sr[good]
min_max_sr = [np.min(sr),np.max(sr)]
fitted_c_200 = np.poly1d(np.polyfit(sr, c_200, 1))(min_max_sr)
print(fitted_c_200)
fitted_c_500 = np.poly1d(np.polyfit(sr, c_500, 1))(min_max_sr)
fitted_c_1000 = np.poly1d(np.polyfit(sr, c_1000, 1))(min_max_sr)
# Fit a linear trend line
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_cadi_1000mas'],'go',label='1000mas')
plt.semilogy(min_max_sr, np.power(10,fitted_c_200),'r--')
plt.semilogy(min_max_sr,np.power(10,fitted_c_500),'b--')
plt.semilogy(min_max_sr,np.power(10,fitted_c_1000),'g--')
plt.xlabel('Measured Strehl')
plt.ylabel('$5\sigma$ post-ADI (classical) contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_contrast_vs_strehl_cadi_with_fit.pdf'))



plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_pca_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_pca_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_pca_1000mas'],'go',label='1000mas')
plt.xlabel('Measured Strehl')
plt.ylabel('$5\sigma$ post-ADI (PCA) contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_contrast_vs_strehl_pca.pdf'))


plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('Measured Strehl')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_strehl.pdf'))


# we add a trend line
c_200 = pd_total['contrast_raw_coro_200mas']
c_500 = pd_total['contrast_raw_coro_500mas']
c_1000 = pd_total['contrast_raw_coro_1000mas']
c_200 = np.log10(c_200[good])
c_500 = np.log10(c_500[good])
c_1000 = np.log10(c_1000[good])
fitted_c_200 = np.poly1d(np.polyfit(sr, c_200, 1))(min_max_sr)
fitted_c_500 = np.poly1d(np.polyfit(sr, c_500, 1))(min_max_sr)
fitted_c_1000 = np.poly1d(np.polyfit(sr, c_1000, 1))(min_max_sr)

plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['measured_strehl'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.semilogy(min_max_sr, np.power(10,fitted_c_200),'r--')
plt.semilogy(min_max_sr,np.power(10,fitted_c_500),'b--')
plt.semilogy(min_max_sr,np.power(10,fitted_c_1000),'g--')
plt.xlabel('Measured Strehl')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_strehl_with_fit.pdf'))


plt.semilogy(pd_total['tau0_O_mean']*1000,pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['tau0_O_mean']*1000,pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['tau0_O_mean']*1000,pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('RTC $\\tau_0$ (ms)')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_RTC_Strehl.pdf'))

plt.semilogy(pd_total['old_dimm_seeing_O_mean'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['old_dimm_seeing_O_mean'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['old_dimm_seeing_O_mean'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('old DIMM seeing (arcsec)')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.axis([0,2.5,3e-5,1e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_old_dimm_seeing.pdf'))

plt.semilogy(pd_total['seeing_O_mean'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['seeing_O_mean'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['seeing_O_mean'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('RTC seeing (arcsec)')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.axis([0,1.,3e-5,1e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_RTC_seeing.pdf'))

plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('RTC turbulent wind speed (m/s)')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
plt.axis([0,25,3e-5,5e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_RTC_windspeed.pdf'))

plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_cadi_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_cadi_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['wind_speed_O_mean'],pd_total['contrast_cadi_1000mas'],'go',label='1000mas')
plt.xlabel('RTC turbulent wind speed (m/s)')
plt.ylabel('$5\sigma$ post-ADI (classical) contrast')
plt.axis([0,25,1e-6,1e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_cadi_contrast_vs_RTC_windspeed.pdf'))

plt.semilogy(pd_total['simbad_FLUX_V'],pd_total['contrast_raw_coro_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['simbad_FLUX_V'],pd_total['contrast_raw_coro_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['simbad_FLUX_V'],pd_total['contrast_raw_coro_1000mas'],'go',label='1000mas')
plt.xlabel('V magnitude')
plt.ylabel('$5\sigma$ raw coronagraphic contrast')
#plt.axis([0,25,3e-5,5e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_raw_coro_contrast_vs_V_mag.pdf'))


plt.semilogy(pd_total['old_dimm_seeing_O_rms'],pd_total['contrast_cadi_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['old_dimm_seeing_O_rms'],pd_total['contrast_cadi_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['old_dimm_seeing_O_rms'],pd_total['contrast_cadi_1000mas'],'go',label='1000mas')
plt.xlabel('RMS seeing (arcsec)')
plt.ylabel('$5\sigma$ post-ADI (classical) contrast')
plt.axis([0,0.1,1e-6,1e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_cadi_contrast_vs_RMS_old_dimm_seeing.pdf'))


plt.semilogy(pd_total['tau0_O_rms'],pd_total['contrast_cadi_200mas'],'ro',label='200mas')
plt.semilogy(pd_total['tau0_O_rms'],pd_total['contrast_cadi_500mas'],'bo',label='500mas')
plt.semilogy(pd_total['tau0_O_rms'],pd_total['contrast_cadi_1000mas'],'go',label='1000mas')
plt.xlabel('RMS  $\\tau_0$ (ms)')
plt.ylabel('$5\sigma$ post-ADI (classical) contrast')
#plt.axis([0,0.1,1e-6,1e-2])
plt.grid()
plt.legend()
plt.savefig(os.path.join(path,'SHARDDS_comparison_cadi_contrast_vs_RMS_tau0_seeing.pdf'))



