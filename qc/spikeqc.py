#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spikeqc.py
Created on Mon May  1 16:40:40 2023

@author: sapje1
"""
import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os

def splash():
    print('\nspikeqc.pc')
    print('   Rudamentory spike quality control check\n')
    print('   USAGE:        python path_to_spikeqc/spikeqc.py SPIKEQC_NIFTI\n')
    print('   ARGUMENTS:    SPIKEQC_NIFTI - spike test nifti')
    return

spike_fname = '/home/sapje1/scratch_sapje1/2023/230430_qc_3TEgradreplacement/spike_check/nifti/DICOM_EPIspike_head_20230428163743_8.nii'
spike_fname = '/home/sapje1/data_sapje1/QC/Connectom/32chcoil.fault/QA3TM_3TW32chcoil_170424/20170524_091834GloverGSQAPQA3TM3TW32chcoils003a001.nii.gz'
spike_fname = '/home/sapje1/data_sapje1/QC/Connectom/32chcoil.fault/QA3TM_3TW32chcoil_170424/20170524_091834GloverGSQAPQA3TM3TW32chcoils012a001.nii.gz'


spike = mriqc.MultiVolQc(spike_fname, False, False)
spike.timeseries(mask=None, plot=True,savepng=False)
print(spike.vol_mean.shape)