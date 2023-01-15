import nibabel as nib
import numpy as np
import os

class fmriqc:
    def __init__(self,path,name):
        self.in_nii_file = name
        self.nii_path = path
        # for testing
        self.nii_path = '/home/sapje1/scratch_sapje1/fmriqc/251_alspac/'
        self.in_nii_file = 'alspacfmri.nii'

    def nii_load(self):
        self.nii_img = nib.load(os.path.join(self.nii_path, self.in_nii_file))
        # this will load as a proxy image
        self.nii_data = self.nii_img.get_fdata()

    def calc_sfnr(self):
        '''
        calc_sfnr():  calculate timeseries mean, standard deviation and signal to
                    fluctuation noise sfnr (Glover)
        '''
        # use same affine transform throughout
        self.affine = self.nii_img.affine
        vol_mean = np.mean(self.nii_data,3)
        # create nifti, using same affine transform as original
        self.nii_mean = nib.Nifti1Image(vol_mean, self.affine)
        vol_stdev = np.std(self.nii_data,3)
        self.nii_stdev = nib.Nifti1Image(vol_stdev, self.affine)
        vol_sfnr = vol_mean / vol_stdev
        self.nii_sfnr = nib.Nifti1Image(vol_sfnr, self.affine)

        nib.save(self.nii_mean, os.path.join(self.nii_path, 'fmriqc_mean.nii'))
        nib.save(self.nii_stdev, os.path.join(self.nii_path, 'fmriqc_stdev.nii'))
        nib.save(self.nii_sfnr, os.path.join(self.nii_path, 'fmriqc_sfnr.nii'))
