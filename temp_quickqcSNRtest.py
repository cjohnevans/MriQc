# -*- coding: utf-8 -*-
"""
This is a temporary script file for testing quickQC NEMA SNR calcs in mriqc
"""
import os, sys
sys.path.append('/home/sapje1/code/MriQc')
import mriqc

os.chdir('/home/sapje1/scratch_sapje1/2024/240220_quickQC/24_02_22-11_12_42-STD-1_3_12_2_1107_5_2_0_19950/nifti')
connectom_files = os.listdir()
f_part1 = '24_02_22-11_12_42-STD-1_3_12_2_1107_5_2_0_19950_'
f_part2 = '_quick_SNR_gre3D.nii'
snr_conn2 = []
snrbgd_conn2 = []
ghost_conn2 = []
for ii in range(2,3):
    f1 = f_part1 + str(ii) + f_part2
    snrcalc = mriqc.BasicQc(write_report=True, \
                            report_path='/home/sapje1/scratch_sapje1/2024/240220_quickQC/report')
    snrcalc.snr_nema_multivol(f1)
    snr_conn2.append(snrcalc.snr)
    snrbgd_conn2.append(snrcalc.snr_bgd)
    ghost_conn2.append(snrcalc.ghost)
    
print(snrcalc.nii_file)

