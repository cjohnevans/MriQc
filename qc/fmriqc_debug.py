#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:40:11 2023

@author: sapje1
"""
import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os

filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi_nifti/20230420121028_QC_MB_GRE_EPI_FA10_Tx200V.nii'

#filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi/short.nii.gz'
data_path = 'C:\\Users\\sapje1\\data\\23_06_16-15_01_26-DST-1_3_12_2_1107_5_2_34_18984'
f = '23_06_16-15_01_26-DST-1_3_12_2_1107_5_2_34_18984_35_Oddity_1p2mm_noMB_TR2_30sl.nii'
filename = os.path.join(data_path,f)
in_vivo = True
fmri = mriqc.FmriQc(filename, in_vivo, run_report=True)

