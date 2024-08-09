#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spikeqc.py
Created on Mon May  1 16:40:40 2023

@author: sapje1
"""
import sys
sys.path.append('/home/sapje1/code/MriQc')
import mriqc

def splash():
    print('\nspikeqc.pc')
    print('   Rudamentory spike quality control check\n')
    print('   USAGE:        python path_to_spikeqc/spikeqc.py SPIKEQC_NIFTI\n')
    print('   ARGUMENTS:    SPIKEQC_NIFTI - spike test nifti')
    return

#spike_fname = '/home/sapje1/scratch_sapje1/2023/230502_qc_spikecodetest/spike_check_3TE/nifti/DICOM_EPIspike_head_20230428163743_8.nii'
#spike_fname = '/home/sapje1/scratch_sapje1/2023/230502_qc_spikecodetest/connectom_bad_spikes/nifti/17_05_24-08_28_57-DST-1_3_12_2_1107_5_2_0_19950_7_EPIspike_head.nii'
spike_fname = '/cubric/collab/108_QA/QA7T/nifti/XNAT_E12622/23_09_28-13_27_06-DST-1_3_12_2_1107_5_2_34_18984_4_spike_readHF.nii'

spike = mriqc.SpikeQc(spike_fname, False, False)
spike.spike_check()
