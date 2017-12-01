#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:34:32 2017
Script that creates the PSF of all files.
@author: jmilli
"""

import os
import numpy as np
import irdisDataHandler as ir
import pdb
import vip
#ds9=vip.fits.vipDS9()
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy import units as u
import pandas as pd
import sphere_utilities as sph
from contrast_utilities import contrast_curve_from_throughput
import image_tools as im

local = True
#local = False

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


target_with_error = []
exception = []

"""
index 5 (HD164249A): The PSF is on the edge of the 200x200 field.
index 11 (HD135599): no O frames the function scaling_factor_psf fails. Marked as bad.
index 12 (FomC): bad conditions.
HD60491: one PSF was removed since the centering failed, but the frame is good
HD38206_2nd_epoch: raw data initially not there, weird...sine the wiki shows the reduced image
HD218340 without PSF analyzed -> no contrast available
HD25457: one PSF with smaller counts rejected.
HD82943: strong low wind effect !

HD14082B_2nd_epoch was not reduced ! Search for companion now
HD22179_2nd_epoch contains only the PSF. Investigate why. I classified it in bad

"""

channels=['left','right']
size_psf = 200
for id_target,target_name in enumerate(targets):
    if target_name!='HD105':
        continue
    pathTarget = os.path.join(pathRoot,target_name)
    pathRaw = os.path.join(pathTarget,'raw')
    pathOut = os.path.join(pathTarget,'pipeline')        
   
#%%
    size = 799
    distarr = im.distance_array((size,size))
    mask = np.logical_and(distarr>100,distarr<350)
    for ichannel,channel in enumerate(channels):
        cube = fits.getdata(os.path.join(pathOut,'HD105_{0:d}x{0:d}_{1:s}_O.fits'.format(size,channel)))
        header = fits.getheader(os.path.join(pathOut,'HD105_{0:d}x{0:d}_{1:s}_O.fits'.format(size,channel)))
        parang = fits.getdata(os.path.join(pathOut,'HD105_parang_O.fits'))
        derot_angles = fits.getdata(os.path.join(pathOut,'HD105_derotation_angles_O.fits'))
    
        nimg = cube.shape[0]
        score = np.ndarray((nimg))
        for i in range(nimg):
            img = cube[i,:,:]*mask
            score[i] = 1./np.sum(img)

        rebin = nimg//100
        binned_indices = im.make_binned_indices(score,rebin,plot=True)
        binned_parang = im.bin_array(parang,binned_indices)
        binned_derot_angles = im.bin_array(derot_angles,binned_indices)
        binned_cube = im.bin_cube(cube,binned_indices)
        nb_binned_frames = len(binned_parang)
        fits.writeto(os.path.join(pathOut,'HD105_{0:d}x{0:d}_{1:s}_O_rebin{2:d}_{3:d}frames.fits'.format(size,channel,rebin,nb_binned_frames)),binned_cube,header=header,clobber=True,output_verify='ignore')
        fits.writeto(os.path.join(pathOut,'HD105_parang_O_rebin{0:d}_{1:d}frames.fits'.format(rebin,nb_binned_frames)),binned_parang,header=header,clobber=True,output_verify='ignore')
        fits.writeto(os.path.join(pathOut,'HD105_derotation_angles_O_rebin{0:d}_{1:d}frames.fits'.format(rebin,nb_binned_frames)),binned_derot_angles,header=header,clobber=True,output_verify='ignore')


        cadi = vip.madi.adi(binned_cube,binned_derot_angles, verbose=False)
        fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi_rebin{3:d}_{4:d}frames.fits'.format(target_name,size,channel,rebin,nb_binned_frames)),cadi,header=header,clobber=True,output_verify='ignore')
        for ncomp in np.arange(2,12,4):
            pca_ncomp = vip.pca.pca(binned_cube,binned_derot_angles, ncomp=ncomp, verbose=True,mask_center_px=8)
            fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_pca_klip_{3:03d}_rebin{4:d}_{5:d}frames.fits'.format(target_name,size,channel,ncomp,rebin,nb_binned_frames)),pca_ncomp,header=header,clobber=True,output_verify='ignore')
