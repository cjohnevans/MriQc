#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:40:11 2023

@author: sapje1
"""
import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc

filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi/20230213_094826QCMBGREEPITx200Vs005a001.nii.gz'
filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi/20230418_141050QCMBGREEPITx200Vs005a001.nii.gz'
filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi_nifti/20230420121028_QC_MB_GRE_EPI_FA10_Tx200V.nii'
#filename = '/home/sapje1/scratch_sapje1/QC/new_fmriqc/mbgreepi/short.nii.gz'
in_vivo = False
fmri = mriqc.FmriQc(filename, in_vivo, run_report=True)
