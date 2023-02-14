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
else:
    nii_path = '/home/sapje1/scratch_sapje1/fmriqc/251_alspac/'
    nii_file = 'alspacfmri.nii'

filename = os.path.join(nii_path, nii_file)

fmri = mriqc.FmriQc(filename, in_vivo=False, run_report=True)


