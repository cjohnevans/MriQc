# -*- coding: utf-8 -*-
"""
This is a temporary script file for testing quickQC NEMA SNR calcs in mriqc
"""
import os, sys
sys.path.append('/home/sapje1/code/MriQc')
import mriqc

path_quick = '/home/sapje1/scratch_sapje1/2024/240220_quickQC/24_02_22-11_12_42-STD-1_3_12_2_1107_5_2_0_19950/nifti'
ff_quick = '24_02_22-11_12_42-STD-1_3_12_2_1107_5_2_0_19950_2_quick_SNR_gre3D.nii'

path_orig= '/home/sapje1/scratch_sapje1/2023/231121_qc_3tw_snr/nema_tse/other_echoes'
ff_orig = '23_09_12-13_23_32-DST-1_3_12_2_1107_5_2_43_66073_4_NEMA_SNR_Uniformity_T2_e2.nii'

qctest = mriqc.BasicQc()
qctest.snr_nema_multivol(os.path.join(path_quick,ff_quick))

