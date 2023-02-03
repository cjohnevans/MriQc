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


fmri = mriqc.fmriqc(nii_path,nii_file,True)
#fmri.calc_sfnr()
fmri.create_report()
#mriqc.plot_histogram(fmri.vol_data,False)
#mriqc.ortho_view(fmri.vol_mean, 'Mean Volume')
#plt.imshow(np.reshape(fmri.vol_mean[:,:,20:41:10].T, (96*3, 96)).T)
#plt.colorbar()
