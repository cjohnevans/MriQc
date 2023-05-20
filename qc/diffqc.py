#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing script for DiffusionQc methods

Created on Sat May 20 06:59:52 2023

@author: sapje1
"""

import sys
import os
sys.path.append('/home/sapje1/python_mrobjects/qc')
import mriqc


prisma_dir = '/home/sapje1/scratch_sapje1/2023/230518_qc_diffusion/prisma/'
prisma_f = '23_05_19-15_52_36-DST-1_3_12_2_1107_5_2_43_66075_3_qcdiff_linear_g_b4000_estG70'
nii_f = os.path.join(prisma_dir,prisma_f+'.nii')

prismaqc = mriqc.MultiVolDiffusion(nii_f)
print(prismaqc.bvals)