#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:40:11 2023

@author: sapje1
"""
import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc

filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi_nifti/20230420121028_QC_MB_GRE_EPI_FA10_Tx200V.nii'
filename = '/home/sapje1/scratch_sapje1/QC/mriqc_7t/nifti/23_05_03-09_47_10-DST-1_3_12_2_1107_5_2_34_18984_6_QC_MB_GRE_EPI_FA15_Tx220V.nii'

#fmri = mriqc.FmriQc(filename)
fmri = mriqc.FmriQc(filename, in_vivo=False, run_report=True)
