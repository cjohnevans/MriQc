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
import matplotlib.pyplot as plt

# wl026
in_dir = '/home/sapje1/scratch_sapje1/2023/230518_qc_diffusion/prisma/'

# capella
in_dir = '/home/john/data/230518_qc_diffusion/prisma'
in_f = '23_05_19-15_52_36-DST-1_3_12_2_1107_5_2_43_66075_3_qcdiff_linear_g_b4000_estG70'

# win10 laptop
in_dir="C:/Users/sapje1/data/diffqc/23_05_18-15_28_59-DST-1_3_12_2_1107_5_2_0_19950"
in_f="23_05_18-15_28_59-DST-1_3_12_2_1107_5_2_0_19950_5_qadiff_linear_g_gmax291"

nii_f = os.path.join(in_dir,in_f+'.nii')

diffqc = mriqc.DiffusionQc(nii_f, nax=3,namp=5, nrep=5,nstart=5,nend=5,b0_in_loop=False)
diffqc.timeseries(plot=True,savepng=False)
plt.show()
diffqc.prep_axis_amp_rep()
diffqc.x.g_uniformity()
