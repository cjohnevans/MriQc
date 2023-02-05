import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os

if os.name == 'nt':
    nii_path = 'C:\\Users\\sapje1\\data\\RoutineQA_examples'
    nii_file = '20211111_113225GloverGSQAPs003a001.nii'
else:
    nii_path = '/home/sapje1/scratch_sapje1/fmriqc/251_alspac/'
    nii_file = 'alspacfmri.nii'


fmri = mriqc.phantomfmriqc(nii_path,nii_file,True)
fmri.voi((10,20,20))
fmri.create_report()
fmri.calc_sfnr(fmri.mask)
