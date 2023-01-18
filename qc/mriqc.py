import nibabel as nib
import numpy as np
import os
import matplotlib.pyplot as plt

# mriqc is generic for checking multi-volume MR data (e.g. fmri, qmt)
class mriqc:
    '''
    fmriqc class: data and methods for dealing with fmri quality control data
    '''
    def __init__(self,path,name, is_invivo):
        self.in_nii_file = name
        self.nii_path = path
        self.invivo = False # default to phantom
        self.sfnr = 0

    def nii_load(self):
        '''
        nii_load()

        Load nifti file using nibabel load.
        Populate nii_img (nibabel) and vol_data (numpy array)
        
        '''        
        self.nii_img = nib.load(os.path.join(self.nii_path, self.in_nii_file))
        # this will load as a proxy image
        self.vol_data = self.nii_img.get_fdata()       
        self.affine = self.nii_img.affine
        if self.vol_data.ndim > 3:
            self.is_multi_volume = True
            self.n_vols = self.vol_data.shape[-1]
        else:
            self.is_multi_volume = False
            self.n_vols = 1

# methods specific to fmri
class fmriqc(mriqc):
    '''
    methods specific to fmri (inherited from mriqc)
    '''
    is_fmri = True

    def calc_sfnr(self):
        '''
        fmriqc.calc_sfnr():  calculate timeseries mean, standard deviation and signal to
                    fluctuation noise sfnr (Glover)

        Populates:
        vol_mean = mean signal across timeseries
        vol_stdev = standard deviation of signal across timepoints
        vol_sfnr = signal to fluctuation noise sfnr (Glover)

        Outputs:
        fmriqc_mean.nii, fmriqc_stdev.nii, fmriqc_sfnr.nii
        
        '''
        # calculate mean, stdev and sfnr
        self.vol_mean = np.mean(self.vol_data,3)
        self.vol_stdev = np.std(self.vol_data,3)
        # mask at 25% of peak voxel intensity of mean image
        self.mask = threshold_vol(self.vol_mean, True, 0.25)
        self.vol_mean = self.vol_mean * self.mask
        self.vol_stdev = self.vol_stdev * self.mask

        # deal with nans.  
        self.vol_sfnr = np.divide(self.vol_mean, self.vol_stdev, \
                             out=np.zeros_like(self.vol_mean), \
                             where=self.vol_stdev!=0)

        ortho_view(self.vol_sfnr, 'SFNR')

        # create nifti, using same affine transform as original
        nii_mean = nib.Nifti1Image(self.vol_mean, self.affine)
        nii_stdev = nib.Nifti1Image(self.vol_stdev, self.affine)
        nii_sfnr = nib.Nifti1Image(self.vol_sfnr, self.affine)
        nib.save(nii_mean, os.path.join(self.nii_path, 'fmriqc_mean.nii'))
        nib.save(nii_stdev, os.path.join(self.nii_path, 'fmriqc_stdev.nii'))
        nib.save(nii_sfnr, os.path.join(self.nii_path, 'fmriqc_sfnr.nii'))


def threshold_vol(vol, by_fraction, threshold):
    '''
    fmriqc.threshold_vol(
       by_fraction
       threshold
    )

    Params:
    by_fraction:  bool.
       True = threshold_val is a fraction of max value in image
       False = threshold_val is absolute pixel intensity
    threshold: float/int
       zero all values below threshold

    Returns:
        np array with mask of values above threshold
    '''

    max_pixel = np.amax(np.amax(np.amax(vol)))
    # mask has to be a copy of vol, to prevent the original volume being overwritten by the mask
    mask = np.array(vol)
    pixel_threshold = max_pixel * threshold
    mask[mask<pixel_threshold] = 0
    mask[mask>=pixel_threshold] = 1
    return mask

def ortho_view(vol, title):
    '''
    mriqc.ortho_view(vol)
s
    Params:
    vol: 3D numpy array to display.  Can't be multi_volume
    
    Show middle slice in three orthogonal views
    '''

    if vol.ndim != 3:
        print('vol has ' + str(vol.ndim) + ' volumes')
        print('Must be a 3D numpy array')
        return

    mid_slice = [int(np.floor(dim_len/2)) for dim_len in vol.shape]
    orth = []
    
    orth.append(vol[:,:,mid_slice[2]])
    orth.append(vol[:,mid_slice[1],:])
    orth.append(vol[mid_slice[0],:,:])
    fig, ax = plt.subplots(1,3)
    for view in [0,1,2]:
        ax[view].imshow(orth[view])
        ax[view].set_xticks([])
        ax[view].set_yticks([])
    ax[1].set_title(title)
 
    return mid_slice

def plot_histogram(vol, save_png):
    '''
    plot_histogram(): 

    Plots a histogram of all pixel values, across all values in 3D or 4D
    numpy array.
    Saves histogram in a png file
    '''
    hist = np.histogram(vol,100)
    bins = hist[1][1:]  # bins includes both ends, take the higher value
    vals = hist[0]
    fig,ax = plt.subplots(1)
    ax.plot(bins, vals)
    ax.set_xlabel('Pixel Value')
    ax.set_ylabel('Number of pixels')
    ax.set_title('Image histogram')
    if save_png:
        fig.savefig(os.path.join(self.nii_path, 'pixel_histogram.png'))
