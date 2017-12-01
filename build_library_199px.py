#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 18:41:04 2017

@author: jmilli
"""

import os
import numpy as np
import pdb
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
import sys
import library as lib

local=False
if local==False:
    library = lib.Library(local=False,channel='left',size=199,pathRDI='RDI_python_dev')
else:
    library = lib.Library(local=True,channel='left',size=199,pathRDI='RDI_python')

#    library.build_covariance_matrix()
#    library.save_covariance_matrix()
library.load_covariance_matrix()
target = 'HD114082'
score = library.analyze_correlation(target,highestest_rank_to_test=50,save=True)
library.build_library(np.where(score>10)[0],filename='library_'+target+'.fits')


