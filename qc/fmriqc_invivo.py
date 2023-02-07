# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 22:02:04 2023

@author: sapje1
"""

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os
import numpy as np
import matplotlib.pyplot as plt

if os.name == 'nt':
    nii_path = 'C:\\Users\\sapje1\\data\\RoutineQA_examples'
    nii_file = 'alspacfmri.nii'
else:
    nii_path = '/home/sapje1/scratch_sapje1/fmriqc/251_alspac/'
    nii_file = 'alspacfmri.nii'    
    
fmri = mriqc.fmriqc(nii_path,nii_file,in_vivo=True)

fmri.create_report()
