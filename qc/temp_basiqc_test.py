# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:07:12 2023

@author: sapje1
"""

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os


#3TW
path = 'C:\\Users\\sapje1\\data\\basicqc\\23_07_24-13_11_02-DST-1_3_12_2_1107_5_2_43_66073//'
f1='scans_NEMA_SNR_1_20230724131352_2_e22.nii'
f2='scans_NEMA_SNR_2_20230724131352_3_e22.nii'
#HF
path = 'C:\\Users\\sapje1\\data\\basicqc\\23_07_24-15_16_01_HG21070042//scans'
f1='scans_T2_(AXI,_Fast)_20230724141652_2.nii'
f2='scans_T2_(AXI,_Fast)_20230724141652_3.nii'
# new HF - K2285 phantom (Siemens cylinder)
path = 'C:\\Users\\sapje1\\data\\basicqc\\hyperfine\\23_07_27-16_03_38_HG21070042'
f1 = '23_07_27-16_03_38_HG21070042_2_T2_(AXI,_Fast)-NOT_FOR_DIAGNOSTIC_USE.nii'
f2 = '23_07_27-16_03_38_HG21070042_3_T2_(AXI,_Fast)-NOT_FOR_DIAGNOSTIC_USE.nii'
f1 = '23_07_27-16_03_38_HG21070042_6_T2_(AXI)-NOT_FOR_DIAGNOSTIC_USE.nii'
f2 = '23_07_27-16_03_38_HG21070042_7_T2_(AXI)-NOT_FOR_DIAGNOSTIC_USE.nii'

# 3TW 230727 - updated protocol
path = 'c:\\Users\\sapje1\\data\\basicqc\\3tw\\23_07_27-09_25_34-DST-1_3_12_2_1107_5_2_43_66073'
f1 = '23_07_27-09_25_34-DST-1_3_12_2_1107_5_2_43_66073_10_NEMA_SNR_Hcoil_imsccorr0p7_3mm_e1.nii'
f2 = '23_07_27-09_25_34-DST-1_3_12_2_1107_5_2_43_66073_11_NEMA_SNR_Hcoil_imsccorr0p7_3mm_e1.nii'
snr = mriqc.BasicQc()
snr.snr_nema(os.path.join(path,f1), os.path.join(path,f2))