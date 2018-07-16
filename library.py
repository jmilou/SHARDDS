#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:34:32 2017
Script that creates the library of frames
@author: jmilli
"""

import os
import numpy as np
#import irdisDataHandler as i
import pdb
#import vip
#ds9=vip.fits.vipDS9()
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
#from astropy import units as u
#import pandas as pd
#import sphere_utilities as sph
from image_tools import distance_array
import matplotlib.gridspec as gridspec 

class Library():
    """This class represents a library
    """

     # class variables

    def __init__(self,local=True,channel='left',size=199,rin=30,rout=60,pathRDI='RDI_python'):
        """
        Class attributes:
            - pathRDI
            - pathRoot
            - pathSparta
            - targets: a list with all targets
            - pathOut_targets: a dictionary with the reduction folder of the targets
            - size: the size of the images we work with
            - channel: the channel we work with (left, right, sum)
            - rin: the inner radius for the correlation
            - rout: the outer radius for the correlation
            - frames_nb: a dictionnary with the number of frames of each target
            - dist_center: an array of shape (size,size) with the distance from the center
            - frames_starta_index:  a dictionnary where keys are target names and values are
                                    the index at which the target starts in the master correlation matrix
            - mask_correlation: an boolean array of shape (size,size) contain 1 for pixels to be used in the correlation
            - total_frames: the total number of frames including all targets
            - correlation_matrix: the large correlation matrix
        """        
        if local:
            self.pathRoot = '/diskb/Data/survey_disk'
            self.pathSparta = '/diskb/Data/sparta'
            self.targets = self.list_targets()
        else:
            self.pathRoot = '/Volumes/SHARDDS_data/survey_disk'
            self.pathSparta = '/Volumes/SHARDDS_data/sparta'
            self.targets = ['HD60491','HD71722','HD218340','HD14082B','HD182681']
        self.pathRDI = os.path.join(self.pathRoot,'../'+pathRDI)

        self.size = size
        self.channel=channel
        self.rin=rin
        self.rout=rout
        self.frames_nb = {}
        self.frames_start_index = {}
        self.frames_start_index_list = np.ndarray((len(self.targets)),dtype=int)       
        self.dist_center = distance_array([self.size,self.size],verbose=False)
#        self.mask_correlation = np.ndarray((self.size,self.size),dtype=int)*0
#        self.mask_correlation[np.logical_and(self.dist_center<rout,self.dist_center>rin)]=1
        self.mask_correlation = np.logical_and(self.dist_center<rout,self.dist_center>rin)
        self.pathOut_targets = {}
        
        total_frames = 0
        for id_target,target_name in enumerate(self.targets):
            pathTarget = os.path.join(self.pathRoot,target_name)
            pathOut = os.path.join(pathTarget,'pipeline')       
            self.pathOut_targets[target_name] = pathOut
            self.frames_start_index[target_name]=total_frames
            self.frames_start_index_list[id_target]=total_frames
            derot_angles = fits.getdata(os.path.join(pathOut,'{0:s}_derotation_angles_O.fits'.format(target_name)))
            nb_frames= len(derot_angles)
            total_frames = total_frames+nb_frames
            self.frames_nb[target_name]=nb_frames
            print('Reading files in {0:s} ({1:d} frames)'.format(target_name,nb_frames))
        self.total_frames = total_frames
        self.correlation_matrix = np.ndarray((total_frames,total_frames),dtype=float)*0.
        
    def list_targets(self):
        """
        """
        targets_tmp = os.listdir(self.pathRoot)
        targets = []        
        for target in targets_tmp:
            if not target.startswith('.f') and os.path.isdir(os.path.join(self.pathRoot,target)) \
            and not target.startswith('DISK') and not target.startswith('bad'):
                targets.append(target)                
        return targets        

    def get_global_index(self,target_name,local_index):
        """
        """
        return self.frames_start_index[target_name]+local_index

    def get_name_and_id_from_index(self,index):
        """
        """
        if index>self.total_frames:
            raise IndexError('Index {0:d} greater than the total number of frames of the library ({1:d})!'.format(index,self.total_frames))
        start_minus_id = index-self.frames_start_index_list
        id_target = len(start_minus_id[start_minus_id>=0])-1
        id_frame = start_minus_id[id_target]
        return self.targets[id_target],id_frame
    
    def get_name_and_id_from_index_list(self,index_list):
        """
        """
        target_list = []
        id_list = []
        for global_index in index_list:
            target,id_frame = self.get_name_and_id_from_index(global_index)
            target_list.append(target)
            id_list.append(id_frame)
        return target_list,id_list

    def build_library(self,index_list,filename=None):
        """
        Input:
            - index_list: a list of indices from the 
            - filename: if not None, saves the library as a fits file with the name
            filename.
        """
        target_list,id_list = self.get_name_and_id_from_index_list(index_list)
        cube = np.ndarray((len(index_list),self.size,self.size))
        unique_targets,target_indices=np.unique(target_list,return_index=False,return_inverse=True,return_counts=False)
        counter = 0
        previous_target_id=-1
        for i,target_index in enumerate(target_indices):
            if previous_target_id != target_index:
                target = target_list[i]
                print('Reading {0:s}'.format(target))
                if self.size==255: # we have to introduce an exception here for the rebinned case.
                    cube_to_load = os.path.join(self.pathOut_targets[target],'{0:s}_1024x1024_rebinned_255x255_{2:s}_O.fits'.format(target,self.size,self.channel))
                else:
                    cube_to_load = os.path.join(self.pathOut_targets[target],'{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target,self.size,self.channel))
                previous_target_id = target_index
                original_cube = fits.getdata(cube_to_load)                
            cube[counter,:,:] = original_cube[id_list[counter],:,:]
            counter=counter+1
        if filename is not None:
            print('Saving the library in {0:s}'.format(os.path.join(self.pathRDI,filename)))
            fits.writeto(os.path.join(self.pathRDI,filename),cube,overwrite=True,output_verify='ignore')
        return cube
            
    def build_covariance_matrix(self):
        """
        """
        for id_target1,target1 in enumerate(self.targets):
            print('Correlating target {0:s} ({1:d}/{2:d})'.format(target1,id_target1+1,len(self.targets)))
            for id_target2 in range(id_target1,len(self.targets)):
                target2  = self.targets[id_target2]
                corr_matrix = self.correlate_targets(target1,target2)
                self.correlation_matrix[self.get_global_index(target1,0):self.get_global_index(target1,self.frames_nb[target1]),\
                                        self.get_global_index(target2,0):self.get_global_index(target2,self.frames_nb[target2])] = \
                                        corr_matrix
                self.correlation_matrix[self.get_global_index(target2,0):self.get_global_index(target2,self.frames_nb[target2]),\
                                        self.get_global_index(target1,0):self.get_global_index(target1,self.frames_nb[target1])] = \
                                        corr_matrix.transpose()
             
    def save_covariance_matrix(self):
        """
        """
        filename = os.path.join(self.pathRDI,'{0:d}x{0:d}_{1:s}_correlation_matrix_rin{2:d}_rout{3:d}.fits'.format(\
                                  self.size,self.channel,self.rin,self.rout))
        print('Saving the correlation matrix in {0:s}'.format(filename))
        fits.writeto(filename,self.correlation_matrix,overwrite=True,output_verify='ignore')
        
    def load_covariance_matrix(self):
        """
        """
        filename = os.path.join(self.pathRDI,'{0:d}x{0:d}_{1:s}_correlation_matrix_rin{2:d}_rout{3:d}.fits'.format(\
                                  self.size,self.channel,self.rin,self.rout))
        print('Loading the correlation matrix from {0:s}'.format(filename))
        self.correlation_matrix = fits.getdata(filename)        
        
    def find_highest_correlated_frames(self,target_name):
        """
        Finds the highest correlated frames for a given target, excluding the 
        target from the library.
        Input:
            - target_name: name of the target
        Output
            - list of indices 
        """
        nframes_target = self.frames_nb[target_name] # number of frames of the target to analyze
        index_target_first_frame = self.get_global_index(target_name,0)
#        index_most_correlated_frames = np.ndarray((nframes_target))       

        target_index_1darray = np.arange(nframes_target) #1D array 
        global_index_1darray = np.arange(self.total_frames) # 1D array
        global_index_subarray, target_index_subarray = np.meshgrid(global_index_1darray,target_index_1darray) #2D array

        index_target_first_frame = self.get_global_index(target_name,0)
        index_target_last_frame  = index_target_first_frame+nframes_target
        sub_correlation_matrix = self.correlation_matrix[index_target_first_frame:index_target_last_frame,:]
        ma_global_index_subarray = np.ma.masked_inside(global_index_subarray,index_target_first_frame,index_target_last_frame-1)

        ma_sub_correlation_matrix = np.ma.masked_array(sub_correlation_matrix,mask=ma_global_index_subarray.mask)        
        correlation_index_max = np.ma.argmax(ma_sub_correlation_matrix,axis=1)        
        return correlation_index_max        

    def correlate_targets(self,target1,target2):
        """
        Computes the correlation of the frames of 2 targets
        """
        target_list = [target1,target2]
        for t in target_list:
            if t not in self.targets:
                print('Target {0:s} no in the target list'.format(t))
                return
        if self.size!=255:
            cube1 = fits.getdata(os.path.join(self.pathOut_targets[target1],'{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target1,self.size,self.channel)))
            cube2 = fits.getdata(os.path.join(self.pathOut_targets[target2],'{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target2,self.size,self.channel)))
        else: # we are in full frame mode here
            cube1 = fits.getdata(os.path.join(self.pathOut_targets[target1],'{0:s}_1024x1024_rebinned_255x255_{1:s}_O.fits'.format(target1,self.channel)))
            cube2 = fits.getdata(os.path.join(self.pathOut_targets[target2],'{0:s}_1024x1024_rebinned_255x255_{1:s}_O.fits'.format(target2,self.channel)))            
        nframes1 = cube1.shape[0]
        nframes2 = cube2.shape[0]
        correlation_matrix = np.ndarray((nframes1,nframes2))
        for iframe1 in range(nframes1):
            for iframe2 in range(nframes2):
#                img1 = cube1[iframe1,:,:]*self.mask_correlation
#                img2 = cube2[iframe2,:,:]*self.mask_correlation
                img1 = cube1[iframe1,:,:][self.mask_correlation]
                img1 -= np.mean(img1)
                img2 = cube2[iframe2,:,:][self.mask_correlation]
                img2 -= np.mean(img2)
                corr_coeff = np.nansum(img1*img2) / ( np.sqrt(np.nansum(img1**2)) * np.sqrt(np.nansum(img2**2)) )
                correlation_matrix[iframe1,iframe2] = corr_coeff
        return correlation_matrix                                
        
    def analyze_correlation(self,target_name,highest_rank_to_test=10,save=True):
        """
        Analyze the correlation of the frames from the library with the frames of
        a specific target. Plots the results in pdf and saves the analysis in a fits file
        Input:
            - target_name: the name of the target to analyze the correlation (string)
            - highest_rank_to_test: for each frame of the target, the algorithm will 
                                    find the correlation coefficient with all the other 
                                    frames of the library. It will sort by descending order
                                    these correlation coefficient and only keep the Nth
                                    frames of the library most highly correlated, N being
                                    the parameter highest_rank_to_test. By default 10.
            - save : 
        Output:
            - score: an 1D array of size the number of frames in the library containing for 
                    each frame the number of times it appears in the top N most 
                    correlated frames with the target frames. This number is therefore 
                    between 0 (it never appears within the top N) and the number of frames
                    of the target (it is always within the top N), where N is 
                    the parameter highest_rank_to_test
        """
        nframes_target = self.frames_nb[target_name] # number of frames of the target to analyze
        nframes_othertargets = self.total_frames-nframes_target # number of remaining frames
        print('The target {0:s} has {1:d} frames, there are {2:d} remaining frames and {3:d} frames in total'.format(\
              target_name,nframes_target,nframes_othertargets,self.total_frames))
        target_index_1darray = np.arange(nframes_target) #1D array 
        global_index_1darray = np.arange(self.total_frames) # 1D array
        global_index_subarray, target_index_subarray = np.meshgrid(global_index_1darray,target_index_1darray) #2D array
#        if save:
#            fits.writeto(os.path.join(self.pathRDI,'global_index.fits'),global_index_subarray,overwrite=True,output_verify='ignore')
#            fits.writeto(os.path.join(self.pathRDI,'target_index.fits'),target_index_subarray,overwrite=True,output_verify='ignore')
        
        # The sub_correlation_matrix has a size (nframes_target,self.total_frames)
        index_target_first_frame = self.get_global_index(target_name,0)
        index_target_last_frame  = index_target_first_frame+nframes_target
        sub_correlation_matrix = self.correlation_matrix[index_target_first_frame:index_target_last_frame,:]
        ma_global_index_subarray = np.ma.masked_inside(global_index_subarray,index_target_first_frame,index_target_last_frame-1)
#        if save:
#            fits.writeto(os.path.join(self.pathRDI,'global_index_mask.fits'),ma_global_index_subarray.mask*1.,overwrite=True,output_verify='ignore')
        ma_sub_correlation_matrix = np.ma.masked_array(sub_correlation_matrix,mask=ma_global_index_subarray.mask)

        if save:
            filename = os.path.join(self.pathRDI,'{0:d}x{0:d}_{1:s}_covariance_matrix_rin{2:d}_rout{3:d}_{4:s}.fits'.format(\
                                      self.size,self.channel,self.rin,self.rout,target_name))
            print('Saving the global covariance matrix in {0:s}'.format(filename))
            fits.writeto(filename,sub_correlation_matrix,overwrite=True,output_verify='ignore')
        
#        # we extract the auto-correlation matrix
#        autocorrelation_matrix  = sub_correlation_matrix[:,self.get_global_index(target_name,0):self.get_global_index(target_name,nframes_target)]
#
#        # we extract the correlation between the target and the other targets
#        if index_target_first_frame ==  0:
#            covariance_matrice = sub_correlation_matrix[:,index_target_last_frame+1:]
#            global_index_1darray_wo_target = np.arange(index_target_last_frame+1,self.total_frames)
#        elif index_target_last_frame ==  self.total_frames-1:
#            covariance_matrice = sub_correlation_matrix[:,0:index_target_first_frame]
#            global_index_1darray_wo_target = np.arange(index_target_first_frame)
#        else:
#            rest1_correlation_matrix = sub_correlation_matrix[:,0:index_target_first_frame]
#            rest2_correlation_matrix = sub_correlation_matrix[:,index_target_last_frame+1:]
#            covariance_matrice = np.concatenate((rest1_correlation_matrix, rest2_correlation_matrix), axis=1)
#            global_index_1darray_wo_target = np.concatenate((np.arange(index_target_first_frame),np.arange(index_target_last_frame+1,self.total_frames)))
#        if save:
#            filename = os.path.join(self.pathRDI,'{0:d}x{0:d}_{1:s}_correlation_matrix_rin{2:d}_rout{3:d}_{4:s}_crosscorrelation.fits'.format(\
#                                      self.size,self.channel,self.rin,self.rout,target_name))
#            print('Saving the target correlation matrix in {0:s}'.format(filename))
#            fits.writeto(filename,autocorrelation_matrix,overwrite=True,output_verify='ignore')
#            filename = os.path.join(self.pathRDI,'{0:d}x{0:d}_{1:s}_correlation_matrix_rin{2:d}_rout{3:d}_{4:s}_covariance.fits'.format(\
#                                      self.size,self.channel,self.rin,self.rout,target_name))
#            print('Saving the covariance matrix limited to other targets in {0:s}'.format(filename))
#            fits.writeto(filename,covariance_matrice,overwrite=True,output_verify='ignore')
            

        # the correlation rank is a matrix of shape (nframes_target,self.total_frames)
        # containing the indices of the target frames (from 0 to nframes_target)
        # sorted such that correlation_rank[:,0] is the indices of the target frames the 
        # the most correlated with each frames from the remaing targets.
        correlation_rank = np.ma.argsort(ma_sub_correlation_matrix,axis=1,fill_value=0)[:,::-1]
        highest_rank_to_test = np.min([highest_rank_to_test,nframes_othertargets])#nframes_target]) # we can't probe more than the number of target frames!
        rank_array = np.arange(1,highest_rank_to_test+1)
        nb_correlated_frames = np.zeros((highest_rank_to_test),dtype=int)
        cum_nb_correlated_frames = np.zeros((highest_rank_to_test),dtype=int)
        min_correlation = np.zeros((highest_rank_to_test))
        min_min_correlation = np.zeros((highest_rank_to_test))
        mean_correlation = np.zeros((highest_rank_to_test))
        indices_correlated_frames = []
        occurence_correlated_frames = []
        score = np.zeros((self.total_frames),dtype=int)
        for irank,rank in enumerate(rank_array):
            highest_correlated_frames,nb_occurences = np.unique(correlation_rank[:,irank],return_counts=True)
            nb_correlated_frames[irank]= len(highest_correlated_frames)
            indices_correlated_frames.append(highest_correlated_frames)
            occurence_correlated_frames.append(nb_occurences)
            if irank in np.arange(3):
                print('The frames {} are the {}th most highly correlated with occurence {}'.format(highest_correlated_frames,rank,nb_occurences))
            for i_id_frame,id_frame in enumerate(highest_correlated_frames):
                score[id_frame] += nb_occurences[i_id_frame]
            min_correlation[irank] = np.min(ma_sub_correlation_matrix[np.arange(nframes_target),correlation_rank[:,irank]])
            min_min_correlation[irank] = np.min(ma_sub_correlation_matrix[:,score>0])            
            mean_correlation[irank] = np.mean(ma_sub_correlation_matrix[:,score>0])            
            cum_nb_correlated_frames[irank] = np.sum(score>0)

        fig = plt.figure(0, figsize=(10,20))
        gs = gridspec.GridSpec(5,1)#height_ratios=[1,1]
        gs.update(left=0.1, right=0.95, bottom=0.1, top=0.93, wspace=0.2, hspace=0.2)
        
        ax1 = plt.subplot(gs[0,0]) # Area for the first plot
        ax2 = plt.subplot(gs[1,0]) # Area for the second plot
        ax3 = plt.subplot(gs[2,0]) # Area for the colorbar
        ax4 = plt.subplot(gs[3,0]) # Area for the colorbar
        ax5 = plt.subplot(gs[4,0]) # Area for the colorbar
        
        ax1.plot(rank_array,cum_nb_correlated_frames)
        ax1.set_ylabel('Number of correlated images')

        ax2.plot(rank_array,min_correlation)
        ax2.set_ylabel('Minimum correlation (frame specific)')

        ax3.plot(rank_array,min_min_correlation)
        ax3.set_ylabel('Minimum correlation (global)')

        ax4.plot(rank_array,mean_correlation)
        ax4.set_ylabel('Mean correlation of selected frames')

        for ax in [ax1,ax2,ax3,ax4]:
            ax.grid(True)
            ax.set_xlabel('Lowest rank probed')

        ax5.semilogy(score,'ro')
        ax5.semilogy([0,self.total_frames-1],np.ones((2))*nframes_target,':b')
        ax5.set_xlabel('Frame number')
        ax5.set_ylabel('Score')
        ax5.grid(True)
        if save:
            plt.savefig(os.path.join(self.pathRDI,\
                '{0:d}x{0:d}_{1:s}_rin{2:d}_rout{3:d}_{4:s}_score_up_to_{5:d}.pdf'.format(\
                 self.size,self.channel,self.rin,self.rout,target_name,highest_rank_to_test)))

#        target_averaged_sub_correlation_matrix = np.mean(sub_correlation_matrix,axis=0)  
#        target_dispersion_sub_correlation_matrix = np.std(sub_correlation_matrix,axis=0)            
#        mask_isTarget = [x in np.arange(self.get_global_index(target_name,0),self.get_global_index(target_name,nframes_target)) for x in range(self.total_frames)]
#        masked_target_averaged_sub_correlation_matrix = np.ma.masked_array(target_averaged_sub_correlation_matrix, mask_isTarget)
#        masked_target_dispersion_sub_correlation_matrix = np.ma.masked_array(target_dispersion_sub_correlation_matrix, mask_isTarget)
#        argmax_covariance,max_covariance = np.argmax(masked_target_averaged_sub_correlation_matrix),np.max(masked_target_averaged_sub_correlation_matrix)
#        max_target_name,max_id = self.get_name_and_id_from_index(argmax_covariance)
#        argmin_covariance,min_covariance = np.argmin(masked_target_averaged_sub_correlation_matrix),np.min(masked_target_averaged_sub_correlation_matrix)
#        min_target_name,min_id = self.get_name_and_id_from_index(argmin_covariance)
#        mean_covariance = np.mean(masked_target_averaged_sub_correlation_matrix)
#
#        print('The maximum covariance for all frames of {0:s} is {1:3.2f} for target {2:s} frame {3:d}'.format(target_name,max_covariance,max_target_name,max_id))
#        print('The minimum covariance for all frames of {0:s} is {1:3.2f} for target {2:s} frame {3:d}'.format(target_name,min_covariance,min_target_name,min_id))            

        plt.figure(1,figsize=(14.17,7))
        gs = gridspec.GridSpec(2,1)#height_ratios=[1,1]
        gs.update(left=0.1, right=0.95, bottom=0.1, top=0.93, wspace=0.2, hspace=0.2)
        
        ax1 = plt.subplot(gs[0,0]) # Area for the first plot
        ax2 = plt.subplot(gs[1,0]) # Area for the second plot
        ax1.plot(np.ma.mean(ma_sub_correlation_matrix,axis=0))
        ax1.set_ylabel('Covariance coefficient')
        ax2.plot(np.ma.std(ma_sub_correlation_matrix,axis=0))
        ax2.set_ylabel('Dispersion of the covariance coefficient')
        for ax in [ax1,ax2]:
            ax.set_xlabel('Frame number (excluding the target)')
            ax.grid(True)
        if save:
            plt.savefig(os.path.join(self.pathRDI,\
                '{0:d}x{0:d}_{1:s}_rin{2:d}_rout{3:d}_{4:s}_mean_std_covariance.pdf'.format(\
                 self.size,self.channel,self.rin,self.rout,target_name)))
        return score
            
if __name__=='__main__':  
    
    library = Library(local=False,channel='left',size=199,pathRDI='RDI_python_dev')     
#    library = Library(local=True,channel='left',size=199,pathRDI='RDI_python')     

#    library.build_covariance_matrix()
#    library.save_covariance_matrix()
    library.load_covariance_matrix()
#    target = 'HD9672'
#    target = 'HD71722'
    target = 'HD105'

#    score = library.analyze_correlation(target,highest_rank_to_test=50,save=True)
#    library.get_name_and_id_from_index(353)
#    library.get_name_and_id_from_index_list([353,340,83])
#    library.build_library([353,340,83],filename='library_HD182681.fits')
    library.build_library(np.where(score>10)[0],filename='library_{0:s}_{1:d}x{1:d}_{2:s}_O.fits'.format(target,library.size,library.channel)

#    highest_correlated_indices = library.find_highest_correlated_frames(target)
#    cube = library.build_library(highest_correlated_indices)
#    cube_target = fits.getdata(os.path.join(library.pathOut_targets[target],target+'_199x199_left_O.fits'))
#    parang_target = fits.getdata(os.path.join(library.pathOut_targets[target],target+'_derotation_angles_O.fits'))
#    import vip
#    ds9=vip.fits.vipDS9()
#    ds9.display(cube,cube_target[0:32,:,:])
#    import adiUtilities as adi
#    cube_subtracted_rdi = adi.simple_reference_subtraction(cube_target,cube,library.mask_correlation)
#    cube_subtracted_adi = adi.subtractMedian(cube_target,cleanMean=1.,mask=None)
#    ds9.display(cube,cube_subtracted_adi,cube_subtracted_rdi)
#    rdi = adi.derotateCollapse(cube_subtracted_rdi,parang_target,rotoff=0.,cleanMean=1.)
#    adi = adi.derotateCollapse(cube_subtracted_adi,parang_target,rotoff=0.,cleanMean=1.)
    
    ds9.display(rdi,adi)
