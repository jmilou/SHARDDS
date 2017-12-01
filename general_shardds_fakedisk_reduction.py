#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:34:32 2017
Script that creates the PSF of all files.
@author: jmilli
"""

import os,sys
sys.path.append('/diska/home/jmilli/lib_py/scattered_light_tool')
sys.path.append('/diska/home/jmilli/lib_py/adi_pipeline')
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
import grater as g
import adiUtilities as adi
import radial_data as rd

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
#size_science = 799
size_science = 199

## to search by index
#id_target = 8
#target_name = targets[id_target]
#print(target_name)

## to search by name
#name = 'HD69830'
#for idx,ta in enumerate(targets):
#    if ta == name:
#        print('The index of {0:s} is {1:d}'.format(name,idx))
#        id_target = idx
#        target_name = targets[id_target]
#print(target_name)

## in case you do a single channel
#ichannel=0
#channel=channels[ichannel]

distance = 10. #pc
a_px = np.array([25,50,75])
a_arcsec = a_px*0.01225
a_au = a_arcsec*distance
inclination = np.array([0,45,75])
nb_inclination = len(inclination)
nb_a = len(a_px)

map_disc_array = np.ndarray((nb_inclination,nb_a,size_science,size_science))    

model = g.Grater(nx=size_science,ny=size_science,distance=distance)
model.set_pa(0)
for i_a,a_au_value in enumerate(a_au):
    for i_incl,incl in enumerate(inclination):
        model.set_inclination(incl)
        model.set_radial_density(ain=5.,aout=-5.,a=a_au_value,e=0.) 
#        disc = model.scattered_light_map
        map_disc_array[i_incl,i_a,:,:] = model.compute_scattered_light()
#fits.writeto(os.path.join(pathFakeDisk,'test_disk.fits'),map_disc_array,clobber=True,output_verify='ignore')    

disc_original_flux = np.ndarray((nb_inclination,nb_a))

for id_target,target_name in enumerate(targets[2:]):
#for id_target,target_name in enumerate(targets[1:]):
    print('Target: {0:s}'.format(target_name))
    try:
        pathTarget = os.path.join(pathRoot,target_name)
        pathRaw = os.path.join(pathTarget,'raw')
        pathOut = os.path.join(pathTarget,'pipeline')        
        pathFakeDisk = os.path.join(pathTarget,'fakeDisk')
        if not os.path.exists(pathFakeDisk):
            os.mkdir(pathFakeDisk)

        fileNames = 'SPHER.*.fits'
        irdis_data = i.IrdisDataHandler(pathRaw,pathOut,fileNames,name=target_name,coordOrigin='derot')        

        px = irdis_data.pixel_scale
        print('Pixel scale {0:6.4f} arcsec/px'.format(px))
        fwhm_theoretical=irdis_data.theoretical_fwhm[0]
        print('Theoretical FWHM {0:4.2f}px'.format(fwhm_theoretical))
        scaling_factor_psf = irdis_data.get_psf_scaling_factor()[0]
        parang,_,_ = irdis_data.get_parang(frameType='O',save=False)
        delta_parang = np.abs(parang[-1]-parang[0])


        psf_data = pd.read_csv(os.path.join(pathOut,'psf_data.csv'))
        psf_sum = fits.getdata(os.path.join(pathOut,'{0:s}_psf_sum.fits'.format(target_name)))
        psf_cropped = (fits.getdata(os.path.join(pathOut,'{0:s}_psf_left_cropped_normalized.fits'.format(target_name)))+\
                       fits.getdata(os.path.join(pathOut,'{0:s}_psf_right_cropped_normalized.fits'.format(target_name))))/2
        psf_cropped = psf_cropped/np.sum(psf_cropped)

        flux_aper = np.sum(psf_data['flux'])*scaling_factor_psf
        flux_max = np.max(psf_sum)*scaling_factor_psf
#        flux_total = np.sum(psf_sum)*scaling_factor_psf

        cube_sum = fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_left_O.fits'.format(\
                    target_name,size_science)))+\
                    fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_right_O.fits'.format(\
                                             target_name,size_science)))
        derotation_angles = fits.getdata(os.path.join(pathOut,'{0:s}_derotation_angles_O.fits'.format(target_name)))

        badframes_left_filename = os.path.join(pathOut,'badframes_left.txt')
        badframes_right_filename = os.path.join(pathOut,'badframes_right.txt')
        badframes_left = np.asarray([line.rstrip('\n') for line in open(badframes_left_filename)],dtype=int)
        badframes_right = np.asarray([line.rstrip('\n') for line in open(badframes_right_filename)],dtype=int)
        badframes = np.unique(np.append(badframes_left,badframes_right))
        goodframes = [index for index in np.arange((derotation_angles.shape[0])) if index not in badframes]

        cube_good = cube_sum[goodframes,:,:]
        derotation_angles_good = derotation_angles[goodframes]

        panda_cadi_sum_contrast = pd.read_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_{3:s}.csv'.format(target_name,size_science,'sum','cadi')))
        cadi_sum = fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_left_cadi.fits'.format(target_name,size_science))) + \
                    fits.getdata(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_right_cadi.fits'.format(target_name,size_science)))
        cadi_rd = rd.Radial_data(cadi_sum)
        cadi_prof = cadi_rd.std       

        # we iterate of the sma first
        
        for i_a,a_px_value in enumerate(a_px):

            id_distance = np.logical_and(panda_cadi_sum_contrast['distance']>a_px_value-2,panda_cadi_sum_contrast['distance']<a_px_value+2)
            noise_aper = np.median(panda_cadi_sum_contrast['noise'][id_distance]* panda_cadi_sum_contrast['sigma corr'][id_distance])
            noise_ADU = noise_aper/flux_aper*flux_max
            print('The noise from vip is {0:.2e} ADU'.format(noise_ADU))
    
            # Reading the cadi image
            cadi_prof_r0_disk = np.mean(cadi_prof[np.logical_and(cadi_rd.r>a_px_value-2,cadi_rd.r<a_px_value+2)])
            print('The noise in the cADi image is {0:.2e} ADU or {1:.2e} of the total PSF'.format(5*cadi_prof_r0_disk,5*cadi_prof_r0_disk/flux_max))
            
            for i_incl,incl in enumerate(inclination):
                
                map_disc = map_disc_array[i_incl,i_a,:,:]
    
                fits.writeto(os.path.join(pathFakeDisk,'{0:s}_{1:d}x{1:d}_initial_normalized_disk_{2:.0f}px_{3:.0f}deg.fits'.format(\
                    target_name,size_science,a_px_value,incl)),map_disc,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
                        
                cubeFakeDisk = adi.insertFakeDisc(np.zeros_like(cube_good),derotation_angles_good,map_disc,rotoff=0.,psf=psf_cropped,copy=False)        
#                cadiFakeDisk = vip.madi.adi(cubeFakeDisk,derotation_angles_good, verbose=False)
                cubeFakeDiskSubtracted = adi.subtractMedian(cubeFakeDisk,cleanMean=1.,mask=None)
                cadiFakeDisk = adi.derotateCollapse(cubeFakeDiskSubtracted,derotation_angles_good,rotoff=0.,cleanMean=0.5)
                throughput_max = np.nanmax(cadiFakeDisk) / np.nanmax(map_disc)
                print('CADI throughput (at maximum disk brightness): {0:.1f}%'.format(throughput_max*100))
                fits.writeto(os.path.join(pathFakeDisk,'{0:s}_{1:d}x{1:d}_cadi_calib_{2:.0f}px_{3:.0f}deg.fits'.format(target_name,size_science,a_px_value,incl)),cadiFakeDisk,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
        
                fakeDisk = map_disc/throughput_max*noise_ADU
                disc_original_flux[i_incl,i_a] = noise_ADU/throughput_max
                cubeFakeDisk = adi.insertFakeDisc(np.zeros_like(cube_good),derotation_angles_good,fakeDisk,rotoff=0.,psf=psf_cropped,copy=False)        
                cubeWithFakeDisk = cube_good+cubeFakeDisk        
#                cubeWithFakeDisk = adi.insertFakeDisc(cube_good,derotation_angles_good,fakeDisk,rotoff=0.,psf=psf_cropped,copy=False)        
                
                cadi = vip.madi.adi(cubeWithFakeDisk,derotation_angles_good, verbose=False)
                fits.writeto(os.path.join(pathFakeDisk,'{0:s}_{1:d}x{1:d}_cadi_{2:.0f}px_{3:.0f}deg.fits'.format(target_name,size_science,a_px_value,incl)),cadi,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
                pca_ncomp = vip.pca.pca(cubeWithFakeDisk,derotation_angles_good, ncomp=10, verbose=True,mask_center_px=8)
                fits.writeto(os.path.join(pathFakeDisk,'{0:s}_{1:d}x{1:d}_pca_klip_10_{2:.0f}px_{3:.0f}deg.fits'.format(target_name,size_science,a_px_value,incl)),pca_ncomp,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')    

        

    except Exception as e:
        print('A problem occured in {0:s}'.format(target_name))
        print(e)
        target_with_error.append(target_name)
        exception.append(e)

if len(target_with_error)>0:
    failedtargets_str = '  \n'.join(target_with_error)
    txtfile = open(os.path.join(pathRoot,'failedtargets_fakedisk.txt'),'w')
    txtfile.write(failedtargets_str)
    txtfile.close()

    