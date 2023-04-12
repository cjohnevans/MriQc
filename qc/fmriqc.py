import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os
import numpy as np
import matplotlib.pyplot as plt

if os.name == 'nt':
    nii_path = 'C:\\Users\\sapje1\\data\\RoutineQA_examples'
    nii_file = '20211111_113225GloverGSQAPs003a001.nii'
    nii_file = 'alspacfmri.nii'
    nii_file = '20211111_113225WarmingUps002a001.nii'
    nii_file = 'WarmingUp_7T.nii.gz'
else:
    nii_path = '/cubric/collab/108_QA/Images/CUBRIC-3TE/QA3TE/23_02_02-14_10_00-STD-1_3_12_2_1107_5_2_43_66075'
    nii_file = 'Warmingup_2.nii'

filename = os.path.join(nii_path, nii_file)

fmri = mriqc.FmriQc(filename, in_vivo=False, run_report=True)


