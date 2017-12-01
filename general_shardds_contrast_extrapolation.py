#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:34:32 2017
Script that creates the PSF of all files.
@author: jmilli
"""

import os
import numpy as np
import irdisDataHandler as i
import pdb
import vip
#ds9=vip.fits.vipDS9()
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
import sphere_utilities as sph
from contrast_utilities import contrast_curve_from_throughput
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from image_tools import distance_array

#local = True
local = False

if local:
    pathRoot = '/diskb/Data/survey_disk'
    pathSparta = '/diskb/Data/sparta'
else:
    pathRoot = '/Volumes/SHARDDS_data/survey_disk'
    pathSparta = '/Volumes/SHARDDS_data/sparta'


#target_names = []
#psf_profile = {}
#contrast = {}
#date = {}
#atm_median_param = {}
#atm_dispersion_param = {}
#parang_var = {}
#paths_target = {}
delta_parang = {}
asymptote = {}
r0 = {}

targets_tmp = os.listdir(pathRoot)
targets = []
for target in targets_tmp:
    if not target.startswith('.f') and os.path.isdir(os.path.join(pathRoot,target)) \
    and not target.startswith('DISK') and not target.startswith('bad'):
        targets.append(target)

#targets=['HD172555']

target_with_error = []
exception = []

channels=['left','right']
size_psf = 200
size_science_fullframe = 799
size_science = 199

## to search by index
#id_target = 8
#target_name = targets[id_target]
#print(target_name)

# to search by name
#name = 'HD181296'
#for idx,ta in enumerate(targets):
#    if ta == name:
#        print('The index of {0:s} is {1:d}'.format(name,idx))
#        id_target = idx
#        target_name = targets[id_target]
#print(target_name)

### in case you do a single channel
#ichannel=0
#channel=channels[ichannel]


def exponential_function(distance, amplitude, r0):
    return amplitude * (1-np.exp(-distance/r0))

zp = 1024#Zero point of the 2MASS Hband filter in Jy
def convert_mag_to_Jy(mag):
    return zp*np.power(10,-mag/2.5)

for id_target,target_name in enumerate(targets):
#for id_target,target_name in enumerate(targets[1:]):
    print('Target: {0:s}'.format(target_name))    
    try:
        pathTarget = os.path.join(pathRoot,target_name)
        pathRaw = os.path.join(pathTarget,'raw')
        pathOut = os.path.join(pathTarget,'pipeline')        
        fileNames = 'SPHER.*.fits'
        irdis_data = i.IrdisDataHandler(pathRaw,pathOut,fileNames,name=target_name,coordOrigin='derot')        
        px = irdis_data.pixel_scale
        if 'HD82943' in target_name:
            star_flux_Jy = 8.39
        else:        
            mag_H = irdis_data.simbad_dico['simbad_FLUX_H']
            star_flux_Jy = convert_mag_to_Jy(mag_H)
        scaling_factor = irdis_data.get_psf_scaling_factor()
        parang,_,_ = irdis_data.get_parang(frameType='O',save=False)
        delta_parang[target_name] = np.abs(parang[-1]-parang[0])
        
          
#%%

        psf_data = pd.read_csv(os.path.join(pathOut,'psf_data.csv'))
        psf_sum = fits.getdata(os.path.join(pathOut,'{0:s}_psf_sum.fits'.format(target_name)))

        distarr_psf  = distance_array((psf_sum.shape[0],psf_sum.shape[1]),verbose=True)
        sky_inner = 70
        sky_outer = 80
        sky_value=np.median(psf_sum[np.logical_and(distarr_psf<sky_outer,distarr_psf>sky_inner)])
        radius_array = np.arange(3,sky_inner)
        psf_flux_array = [np.nansum(psf_sum[distarr_psf<r]-sky_value) for r in radius_array]
        plt.close()
        plt.plot(radius_array,psf_flux_array)
        plt.xlabel('Radius in px')
        plt.ylabel('Integrated flux in ADU')
        plt.savefig(os.path.join(pathOut,'PSF_integrated_flux.pdf'))
        plt.close()

        ref_flux = np.max(psf_flux_array)*scaling_factor[0]
        
        red_types = ['pca_klip_010','cadi']
        for ired,red_type in enumerate(red_types): 
            cadi_sum = fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_{3:s}.fits'.format(target_name,size_science,\
                                    channels[0],red_type))) + fits.getdata(os.path.join(pathOut,\
                                    '{0:s}_{1:d}x{1:d}_{2:s}_{3:s}.fits'.format(target_name,size_science,channels[1],red_type)))
            fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_{3:s}.fits'.format(target_name,size_science,'sum',red_type)),\
                         cadi_sum,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
            
            panda_contrast_curve_left = pd.read_csv(os.path.join(pathOut,\
                            '{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_{3:s}.csv'.format(target_name,size_science,channels[0],red_type)))
            popt_left, pcov_left = curve_fit(exponential_function, panda_contrast_curve_left['distance'], \
                                             panda_contrast_curve_left['throughput'],bounds=([0.,10.],[1.,2000.]))    
            panda_contrast_curve_right = pd.read_csv(os.path.join(pathOut,\
                            '{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_{3:s}.csv'.format(target_name,size_science,channels[1],red_type)))
            popt_right, pcov_right = curve_fit(exponential_function, panda_contrast_curve_right['distance'], \
                                               panda_contrast_curve_right['throughput'],bounds=([0.,10.],[1.,2000.]))
            start_distance = np.fix(np.min(np.append(panda_contrast_curve_left['distance'],panda_contrast_curve_right['distance'])))+1
            end_distance = np.fix(np.max(np.append(panda_contrast_curve_left['distance'],panda_contrast_curve_right['distance'])))+1
            mean_distance = np.arange(start_distance,end_distance)
            interp_function_left = interp1d(panda_contrast_curve_left['distance'],panda_contrast_curve_left['throughput'],\
                               kind='linear',bounds_error=False,fill_value='extrapolate')            
            throughput_interp_left = interp_function_left(mean_distance) 
            interp_function_right = interp1d(panda_contrast_curve_right['distance'],panda_contrast_curve_right['throughput'],\
                               kind='linear',bounds_error=False,fill_value='extrapolate')           
            throughput_interp_right = interp_function_right(mean_distance) 
            mean_throughput = (throughput_interp_left+throughput_interp_right)/2.
    
            panda_cadi_sum_contrast = contrast_curve_from_throughput(cadi_sum,\
                            np.mean(psf_data['fwhm']), px,\
                            np.sum(psf_data['flux'])*psf_data['scaling'][0],\
                            throughput=(mean_distance,mean_throughput),sigma=5)
            panda_cadi_sum_contrast['sensitivity (Student) [mJy/arcsec^2]'] = panda_cadi_sum_contrast['sensitivity (Student)']*\
                    (np.sum(psf_data['flux'])*psf_data['scaling'][0])/ref_flux*star_flux_Jy*1000/(px**2)
            panda_cadi_sum_contrast.to_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_{3:s}.csv'.format(target_name,size_science,'sum',red_type)))
            
#            plt.close()
#            plt.figure(ired)
#            plt.semilogy(panda_contrast_curve_left['distance'],panda_contrast_curve_left['sensitivity (Student)'], 'r-', label="Left")
#            plt.semilogy(panda_contrast_curve_right['distance'],panda_contrast_curve_right['sensitivity (Student)'], 'b-', label="Right")
#            plt.semilogy(panda_cadi_sum_contrast['distance'],panda_cadi_sum_contrast['sensitivity (Student)'], 'g-', label="Sum")
#            plt.legend(frameon=False)
#            plt.xlabel('Separation in px')
#            plt.ylabel('Contrast')
#            plt.savefig(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_{3:s}.pdf'.format(target_name,size_science,'sum',red_type)))        

#        plt.close()
#        plt.figure(0)                
#        plt.plot(panda_contrast_curve_left['distance'],panda_contrast_curve_left['throughput'], 'ro', label="Measured throughput left")
#        plt.plot(panda_contrast_curve_left['distance'], exponential_function(panda_contrast_curve_left['distance'], *popt_left), 'r-', label="Fitted throughput left")
#        plt.plot(panda_contrast_curve_right['distance'],panda_contrast_curve_right['throughput'], 'bo', label="Measured throughput right")
#        plt.plot(panda_contrast_curve_right['distance'], exponential_function(panda_contrast_curve_right['distance'], *popt_right), 'b-', label="Fitted throughput right")        
#        plt.legend(frameon=False)

        asymptote[target_name] = (popt_left[0]+popt_right[0])/2.
        r0[target_name] = (popt_left[1]+popt_right[1])/2.

        cadi_fullframe_name_left = os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_fullframe,channels[0]))
        cadi_fullframe_name_right = os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_fullframe,channels[1]))

#        if (os.path.isfile(cadi_fullframe_name_left) and os.path.isfile(cadi_fullframe_name_right)) == False:
#            continue
#        else:
        cadi_sum_fullframe = fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_fullframe,\
                                channels[0]))) + fits.getdata(os.path.join(pathOut,\
                                '{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_fullframe,channels[1])))
        fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_fullframe,'sum')),cadi_sum_fullframe,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
        
        distance_extrapolated_fullframe = np.arange(np.max(mean_distance)+1,size_science_fullframe/2.)
        distance_fullframe = np.append(mean_distance,distance_extrapolated_fullframe)
        throughput_fullframe = np.append(mean_throughput,exponential_function(distance_extrapolated_fullframe, *(popt_right+popt_left)/2.))

#        plt.close()
#        plt.semilogx(panda_contrast_curve_left['distance'],panda_contrast_curve_left['throughput'], 'ro', label="Measured throughput left")
##        plt.semilogx(panda_contrast_curve_left['distance'], exponential_function(panda_contrast_curve_left['distance'], *popt_left), 'r-', label="Fitted throughput left")
#        plt.semilogx(panda_contrast_curve_right['distance'],panda_contrast_curve_right['throughput'], 'bo', label="Measured throughput right")
##        plt.semilogx(panda_contrast_curve_right['distance'], exponential_function(panda_contrast_curve_right['distance'], *popt_right), 'b-', label="Fitted throughput right")
#        plt.semilogx(distance_fullframe,throughput_fullframe,label="Extrapolated throughput sum")
#        plt.legend(frameon=False)        
#        plt.xlabel('Separation in px')
#        plt.ylabel('Throughput')
#        plt.savefig(os.path.join(pathOut,'{0:s}_extrapolated_throughput_{1:s}.pdf'.format(target_name,red_type)))        
        
        panda_cadi_sum_fullframe_contrast = contrast_curve_from_throughput(cadi_sum_fullframe,\
                        np.mean(psf_data['fwhm']), px,\
                        np.sum(psf_data['flux'])*psf_data['scaling'][0],\
                        throughput=(distance_fullframe,throughput_fullframe),sigma=5)
        panda_cadi_sum_fullframe_contrast['sensitivity (Student) [mJy/arcsec^2]'] = panda_cadi_sum_fullframe_contrast['sensitivity (Student)']*\
                (np.sum(psf_data['flux'])*psf_data['scaling'][0])/ref_flux*star_flux_Jy*1000/(px**2)
        panda_cadi_sum_fullframe_contrast.to_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_sum_sorted_cadi.csv'.format(target_name,size_science_fullframe)))
        

#%%  rebin by a facto 4 of the images
    
    except Exception as e:
        print('A problem occured in {0:s}'.format(target_name))
        print(e)
        target_with_error.append(target_name)
        exception.append(e)

if len(target_with_error)>0:
    failedtargets_str = '  \n'.join(target_with_error)
    txtfile = open(os.path.join(pathRoot,'failedtargets.txt'),'w')
    txtfile.write(failedtargets_str)
    txtfile.close()


#plt.scatter(np.asarray(delta_parang.values()),np.asarray(asymptote.values()))
#plt.ylabel('Highest throughput')
#plt.xlabel('Amplitude of parallactic angle variation')
#plt.grid()
#plt.savefig('/Users/jmilli/Desktop/ADI_max_throughput.pdf')
#
#r0_resel = np.asarray(r0.values())/4.
#plt.semilogx(np.asarray(delta_parang.values()),r0_resel,'bo')
#plt.ylabel('Separation in resel to achieve 1/e throughput')
#plt.xlabel('Amplitude of parallactic angle variation in $\circ$')
#plt.grid()
#plt.savefig('/Users/jmilli/Desktop/ADI_throughput_separation.pdf')
#
#nb_resel_scanned_at_r0 = r0_resel*np.deg2rad(np.asarray(delta_parang.values()))
#plt.scatter(np.asarray(delta_parang.values()),nb_resel_scanned_at_r0)
#plt.ylabel('Rotation in resels to achieve 1/e throughput')
#plt.xlabel('Amplitude of parallactic angle variation')
#plt.grid()
#plt.savefig('/Users/jmilli/Desktop/ADI_throughput_rotation.pdf')
