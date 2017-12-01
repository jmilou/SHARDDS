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

local = True
#local = False

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
#
#targets_tmp = os.listdir(pathRoot)
#targets = []
#for target in targets_tmp:
#    if not target.startswith('.f') and os.path.isdir(os.path.join(pathRoot,target)) \
#    and not target.startswith('DISK') and not target.startswith('bad'):
#        targets.append(target)
#
#
#target_with_error = []
#exception = []

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
#for id_target,target_name in enumerate(targets[0:1]):

target_name = 'HD14082B_2nd_epoch'

pathTarget = os.path.join(pathRoot,target_name)
pathRaw = os.path.join(pathTarget,'raw')
pathOut = os.path.join(pathTarget,'pipeline')

fileNames = 'SPHER.*.fits'
irdis_data = i.IrdisDataHandler(pathRaw,pathOut,fileNames,name=target_name,coordOrigin='derot')

px = irdis_data.pixel_scale
print('Pixel scale {0:6.4f} arcsec/px'.format(px))
fwhm_theoretical=irdis_data.theoretical_fwhm[0]
print('Theoretical FWHM {0:4.2f}px'.format(fwhm_theoretical))

pathSpartaNight = os.path.join(pathSparta,str((irdis_data.date_start-12*u.hour).datetime.date())) 
irdis_data.analyse_sparta(pathSpartaNight,debug=False,force=False)

irdis_data.writeMetaData()

irdis_data.compute_statistics(frameType='all')
irdis_data.compute_statistics(frameType='O')
irdis_data.compute_statistics(frameType='F')


#%%
# We format the PSF 
#        
dico_psf = {'fwhmx':[],'fwhmy':[],'fwhm':[],'flux':[],\
                'scaling':irdis_data.get_psf_scaling_factor(),\
                'fwhm_diffraction':irdis_data.theoretical_fwhm}
    
for ichannel,channel in enumerate(channels):
    psf_cube = irdis_data.get_psf_frames(camera=channel,size=size_psf)
    #ds9.display(psf_cube)
    psf = np.median(psf_cube,axis=0)
    fits.writeto(os.path.join(pathOut,'{0:s}_psf_{1:s}.fits'.format(target_name,channel)),psf,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
    #ds9.display(psf)
    
    centerx,centery,fwhmx,fwhmy,amp,theta = vip.var.fit_2dgaussian(psf, \
        crop=True, cropsize=9, debug=True,full_output=True)
    fwhm_measured = np.mean([fwhmx,fwhmy])
    print('Measured FWHM (left) {0:4.2f}'.format(fwhm_measured)) #4.04
    dico_psf['fwhmx'].append(fwhmx)
    dico_psf['fwhmy'].append(fwhmy)
    dico_psf['fwhm'].append(fwhm_measured)
    
    nx_psf,ny_psf=psf.shape
    psf_flux = vip.phot.aperture_flux(psf, [ny_psf/2], [nx_psf/2], fwhm_measured, ap_factor=1)
    print('Flux of the {0:s} PSF in 1 FWHM: {1:e}'.format(channel,psf_flux[0]))
    dico_psf['flux'].append(psf_flux[0])
    
    size_cropped_psf=10
    psf_c = psf[nx_psf/2-size_cropped_psf:nx_psf/2+size_cropped_psf+1,ny_psf/2-size_cropped_psf:ny_psf/2+size_cropped_psf+1]
    psf_c[psf_c<0]=0
    psf_c = psf_c / psf_flux[0]
    #ds9.display(psf_c)
    fits.writeto(os.path.join(pathOut,'{0:s}_psf_{1:s}_cropped_normalized.fits'.format(target_name,channel)),\
                 psf_c,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
    
psf_data = pd.DataFrame(dico_psf)
psf_data.to_csv(os.path.join(pathOut,'psf_data.csv'),index=False)
psf_sum = fits.getdata(os.path.join(pathOut,'{0:s}_psf_{1:s}.fits'.format(target_name,channels[0])))+\
           fits.getdata(os.path.join(pathOut,'{0:s}_psf_{1:s}.fits'.format(target_name,channels[1])))
fits.writeto(os.path.join(pathOut,'{0:s}_psf_sum.fits'.format(target_name)),psf_sum,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')

   

#%%
size_science = 199
for ichannel,channel in enumerate(channels):
    cube_center,_,_ = irdis_data.write_master_cube(camera=channel,size = size_science,frameType='C',output=True)
    #fits.writeto(os.path.join(pathOut,'cube_waffle_left_before_centering.fits'),cube_center_left,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
    #ds9.display(cube_center)
    
    center,shift = sph.sph_irdis_centerspot(cube_center, filter_name='B_H', \
        centerguessxy=np.ones(2)*size_science//2,path=pathOut,name='waffle_{0:s}'.format(channel),sigfactor=8) #sigfactor=6  centerguessxy=center_xy[ichannel]
    print(center)
    print(shift)

    shift_mean = np.mean(shift,axis=0)

    cube_center_recentered,_,_ = irdis_data.write_master_cube(camera=channel,\
        centerxy=irdis_data._centerxy[ichannel,:]+shift_mean,size=size_science,frameType='C',output=True)
     #ds9.display(cube_center_recentered)

    cube,parang,derot_angles = irdis_data.write_master_cube(camera=channel,\
        centerxy=irdis_data._centerxy[ichannel,:]+shift_mean,size=size_science,frameType='O',output=True)
    #ds9.display(cube)

    goodframes,badframes = vip.preproc.badframes.cube_detect_badfr_pxstats(cube, \
        mode='annulus', in_radius=23, width=30, top_sigma=2., low_sigma=2.0,\
        window=None,plot=True, verbose=True)   
    badframes_str = '  \n'.join([str(bad) for bad in badframes])
    txtfile = open(os.path.join(pathOut,'badframes_{0:s}.txt'.format(channel)),'w')
    txtfile.write(badframes_str)
    txtfile.close()
    
    cube_good = cube[goodframes,:,:]
    parang_good = parang[goodframes]
    derot_angles_good = derot_angles[goodframes]
    fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_O_sorted.fits'.format(target_name,size_science,channel)),cube_good,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
    fits.writeto(os.path.join(pathOut,'{0:s}_parang_{1:s}_O_sorted.fits'.format(target_name,channel)),parang_good,clobber=True,output_verify='ignore')
    fits.writeto(os.path.join(pathOut,'{0:s}_derotation_angles_{1:s}_O_sorted.fits'.format(target_name,channel)),derot_angles,clobber=True,output_verify='ignore')


    cadi = vip.madi.adi(cube_good,derot_angles_good, verbose=False)
    fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_cadi.fits'.format(target_name,size_science,channel)),cadi,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')
    for ncomp in np.arange(2,12,4):
        pca_ncomp = vip.pca.pca(cube_good,derot_angles_good, ncomp=ncomp, verbose=True,mask_center_px=8)
        fits.writeto(os.path.join(pathOut,'{0:s}_{1:d}x{1:d}_{2:s}_pca_klip_{3:03d}.fits'.format(target_name,size_science,channel,ncomp)),pca_ncomp,header=irdis_data.firstHeader,clobber=True,output_verify='ignore')


    psf_c = fits.getdata(os.path.join(pathOut,'{0:s}_psf_{1:s}_cropped_normalized.fits'.format(target_name,channel)))
    psf_data = pd.read_csv(os.path.join(pathOut,'psf_data.csv'))
    panda_contrast_curve = \
        vip.phot.contrast_curve(cube_good,derot_angles_good, \
        psf_c, psf_data['fwhm'][ichannel],px, \
        psf_data['flux'][ichannel]*psf_data['scaling'][ichannel],\
        vip.madi.adi,sigma=5,nbranch=3, debug=True,\
        plot=True,student=True,verbose=True)
    panda_contrast_curve.to_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_cadi.csv'.format(target_name,size_science,channel)))

    panda_contrast_curve = \
        vip.phot.contrast_curve(cube_good,derot_angles_good, \
        psf_c, psf_data['fwhm'][ichannel],px, \
        psf_data['flux'][ichannel]*psf_data['scaling'][ichannel],\
        vip.pca.pca,sigma=5,nbranch=3, debug=True,\
        plot=True,student=True,verbose=True,mask_center_px=8,ncomp=ncomp)
    panda_contrast_curve.to_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_pca_klip_{3:03d}.csv'.format(target_name,size_science,channel,ncomp)))

    median_cube_good = np.median(cube_good,axis=0)
    panda_raw_contrast = contrast_curve_from_throughput(median_cube_good,\
                    psf_data['fwhm'][ichannel], px,\
                    psf_data['flux'][ichannel]*psf_data['scaling'][ichannel],\
                    throughput=None,sigma=5)
    panda_raw_contrast.to_csv(os.path.join(pathOut,'{0:s}_contrast_curve_{1:d}x{1:d}_{2:s}_sorted_raw_coro.csv'.format(target_name,size_science,channel)))


