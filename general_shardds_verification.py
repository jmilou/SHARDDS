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

#local = True
local = False

if local:
    pathRoot = '/diskb/Data/survey_disk'
    pathSparta = '/diskb/Data/sparta'
else:
    pathRoot = '/Volumes/SHARDDS_data/survey_disk'
    pathSparta = '/Volumes/SHARDDS_data/sparta'


target_names = []
psf_profile = {}
contrast = {}
date = {}
atm_median_param = {}
atm_dispersion_param = {}
parang_var = {}
paths_target = {}

targets_tmp = os.listdir(pathRoot)
targets = []
for target in targets_tmp:
    if not target.startswith('.f') and os.path.isdir(os.path.join(pathRoot,target)) \
    and not target.startswith('DISK') and not target.startswith('bad'):
        targets.append(target)

missing_files = []
target_with_error = []
exception = []

channels=['left','right']
size_psf = 200
size_science_small = 199
size_science_full = 799

## to search by index
#id_target = 8
#target_name = targets[id_target]
#print(target_name)

## to search by name
#name = 'HD105'
#for idx,ta in enumerate(targets):
#    if ta == name:
#        print('The index of {0:s} is {1:d}'.format(name,idx))
#        id_target = idx
#        target_name = targets[id_target]
#print(target_name)
##
## in case you do a single channel
#ichannel=0
#channel=channels[ichannel]

for id_target,target_name in enumerate(targets):
    print('Target: {0:s}'.format(target_name))
    try:
        pathTarget = os.path.join(pathRoot,target_name)
        pathRaw = os.path.join(pathTarget,'raw')
        pathOut = os.path.join(pathTarget,'pipeline')        
    
        names_to_test = []
        for ichannel,channel in enumerate(channels):
            names_to_test.append(os.path.join(pathOut,'{0:s}_psf_{1:s}.fits'.format(target_name,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_psf_{1:s}_cropped_normalized.fits'.format(target_name,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target_name,size_science_full,channel)))            
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_O_sorted.fits'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science_full,channel)))
            for ncomp in np.arange(2,12,4):
                names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_pca_klip_{3:03d}.fits'.format(target_name,size_science_small,channel,ncomp)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_pca_klip_{3:03d}.csv'.format(target_name,size_science_small,channel,ncomp)))
            names_to_test.append(os.path.join(pathOut,'badframes_{0:s}_fullframe.txt'.format(channel))) # not done yet
            names_to_test.append(os.path.join(pathOut,'badframes_{0:s}.txt'.format(channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_cadi.csv'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_raw_coro.csv'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_O_sorted.fits'.format(target_name,size_science_small,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_parang_{1:s}_O_sorted.fits'.format(target_name,channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_derotation_angles_{1:s}_O_sorted.fits'.format(target_name,channel)))   
            derotation_angles_filename = os.path.join(pathOut,'{0:s}_derotation_angles_O.fits'.format(target_name))
            names_to_test.append(derotation_angles_filename)  
            derotation_angles_hdu = fits.open(derotation_angles_filename)
            length_derotation_angle = derotation_angles_hdu[0].header['NAXIS1']
            derotation_angles_hdu.close()
#            print(length_derotation_angle)
            rebinned_cube_filename = os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_{2:s}_O.fits'.format(target_name,1024,channel))
            names_to_test.append(rebinned_cube_filename)
            rebinned_cube_hdu = fits.open(rebinned_cube_filename)
            length_rebinned_cube = rebinned_cube_hdu[0].header['NAXIS3']
            rebinned_cube_hdu.close()
#            print(length_rebinned_cube)
            if length_derotation_angle != length_rebinned_cube:
                print('Problem with {0:s}: rebinned cube and derotation angle have a different size'.format(target_name))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_{2:s}_cadi.fits'.format(target_name,1024,channel)))
            for ncomp in np.arange(2,12,4):
                names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_{2:s}_pca_klip_{3:03d}.fits'.format(target_name,1024,channel,ncomp)))
                names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_{2:s}_pca_klip_{3:03d}_filtered.fits'.format(target_name,1024,channel,ncomp)))
            names_to_test.append( os.path.join(pathOut,'badframes_rebin_{0:s}.txt'.format(channel)))
            names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_{2:s}_cadi_filtered.fits'.format(target_name,1024,channel)))
        names_to_test.append(os.path.join(pathOut,'{0:s}_psf_sum.fits'.format(target_name)))
        names_to_test.append(os.path.join(pathOut,'psf_data.csv'))
#        names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_O.fits'.format(target_name,1024)))# not there yet
#        names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_sum_cadi.fits'.format(target_name,1024))) # not there yet
#        names_to_test.append(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_rebinned_255x255_sum_cadi_filtered.fits'.format(target_name,1024)))

        for names in names_to_test:
            if not os.path.isfile(names): 
                missing_files.append(names)
            
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
else:
    print('No errors')

if len(missing_files)>0:
    missingfiles_str = '  \n'.join(missing_files)
    txtfile = open(os.path.join(pathRoot,'missingfiles.txt'),'w')
    txtfile.write(missingfiles_str)
    txtfile.close()
else:
    print('No missing files')

     