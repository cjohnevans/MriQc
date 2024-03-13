import nibabel as nib
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import datetime as dt
import time  # for dealing with linux file timestamps in FmriQcOverview
import subprocess

class BasicQc:
    '''
    BasicQc: methods for elementary QC calculations
    '''
    def __init__(self, write_report=False, report_path='.'):
        self.write_report = write_report 
        self.report_path = report_path
        if not os.path.exists(self.report_path):
            os.mkdir(self.report_path)
        self.nii_file = '' # set this when loading
        
    # Helper functions first
    def locate_phantom(self,im1):
        """
        locate_phantom(im1)

        Parameters
        ----------
        im1 : image data
            DESCRIPTION.

        Returns
        -------
        (image_centre_x, image_centre_y), (phantom_centre_x, phantom_centre_y).

        """
        # find image centre
        self.im_centre = [int(i/2) for i in im1.shape]
        # find phantom centre, use im1
        # first get the max extent of the phantom,  
        # max collapses along the axis defined by axis               
        col_max = np.ravel(np.max(im1,axis=0)) # collapse 
        row_max = np.ravel(np.max(im1,axis=1))
        
        # then find the edges ...
        col_diff = np.absolute(np.ravel(np.diff(col_max)))
        row_diff = np.absolute(np.ravel(np.diff(row_max)))
        col_edge1 = np.argmax(col_diff[0:self.im_centre[1]])
        col_edge2 = self.im_centre[1]+np.argmax(col_diff[ (self.im_centre[1]+1):] )
        # ... to get the centre
        col_centre = (col_edge1 + col_edge2)/2
        col_width = col_edge2 - col_edge1
        row_edge1 = np.argmax(row_diff[0:self.im_centre[0]])
        row_edge2 = self.im_centre[0]+np.argmax(row_diff[ (self.im_centre[0]+1):] )        
        row_centre = (row_edge1 + row_edge2)/2
        row_width = row_edge2 - row_edge1
        self.ph_centre = (row_centre,col_centre)
        # calc minimum expected radius
        self.ph_radius = np.min([row_width, col_width])/2
    
    def mask_phantom(self, R_mask=None, plot_mask=False):
        '''
        mask_phantom:  create mask of the phantom - i.e. mask out the exterior
        '''
        mask = self.mask_image(region='in', R_mask=R_mask, plot_mask=plot_mask)
        return mask
    
    def mask_empty(self, R_mask=None, plot_mask=False):
        '''
        mask_empty: create mask of empty space around phantom (for noise calcs)
        '''
        mask = self.mask_image(region='out', R_mask=R_mask, plot_mask=plot_mask)
        return mask
        
    def mask_image(self, region='in', R_mask=None, plot_mask=False):
        '''
        mask_image(region='in',
                   R_mask=None, 
                   plot_mask=False)

        Parameters
        ----------
        region : string, required
            accepted values are;
            'in' - to include everything inside R_mask within the mask . 
            'out'- to include everything outside R_mask within the mask . 

            The default is 'in'.
        R_mask : FLOAT, optional
            DESCRIPTION. The default is None.
        plot_mask : BOOL, optional
            Create plot of the mask. The default is False.

        Returns
        -------
        mask : BOOL
            Boolean image mask.

        '''
        row_min = int(-self.ph_centre[0])
        row_max = int(self.im_shape[0] - self.ph_centre[0])
        col_min = int(-self.ph_centre[1])
        col_max = int(self.im_shape[1] - self.ph_centre[1])

        # work out x^2+y^2 <= R^2 using a meshgrid
        # watch meshgrid introduces a rotation
        col_grid, row_grid = np.meshgrid(range(col_min,col_max+1), range(row_min,row_max+1))
        col_grid2 = col_grid[0:self.im1.shape[0], 0:self.im1.shape[1]]
        row_grid2 = row_grid[0:self.im1.shape[0], 0:self.im1.shape[1]]
        col_sqr = np.power(col_grid2,2)
        row_sqr = np.power(row_grid2,2)
        if region == 'out':
            mask = (col_sqr + row_sqr) > np.power(R_mask,2)
        else: # region == 'in'
            mask = (col_sqr + row_sqr) < np.power(R_mask,2)
        
        
        if plot_mask:
            mask_img = mask*np.max(np.max(self.im1))
            fig3 = plt.figure()
            ax = fig3.subplots(1,1)       
            ax.imshow(mask_img+self.im1)
            ax.set_title('mask')
            
        return mask
    
    def mask_corners(self, corner_length=20, plot_mask=False):
        '''
        mask_corners(corner_length, 
                     plot_mask)

        Parameters
        ----------
        corner_length : FLOAT, required
            Length, in pixels, occupied by each corner.
            The default is 20.
        plot_mask : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        mask

        '''
        
        # use the shape of one of the image files - i.e. after converting to 2D
        mask = np.empty(self.im1.shape, dtype='bool')
        mask[:] = False
        n0=self.im1.shape[0]
        n1=self.im1.shape[1]
        mask[:corner_length, :corner_length] = True
        mask[(n0-corner_length):, :corner_length ] = True
        mask[:corner_length, (n1-corner_length):] = True
        mask[(n0-corner_length):,(n1-corner_length):] = True
        # print('corners :', corner_length, n0-corner_length, n1-corner_length)
        
        if plot_mask:
            mask_img = mask
            fig3 = plt.figure()
            ax = fig3.subplots(1,1)       
            ax.imshow(mask_img)
            ax.set_title('mask')
        
        return mask
    
    # calculation functions
    def snr_nema(self, f1, f2):
        """
        snr_nema(nii1,
                 nii2)

        Orientation:  axis 0 is the vertical direction (y) 
                      axis 1 is the horizontal direction (x) 
                      both for the referencing of the array and on imshow
                      with im[axis0,axis1]


        Parameters
        ----------
        f1 : string
            path to first acquisition used in SNR calculation.
        f2 : string
            path to first acquisition used in SNR calculation.

        Returns
        -------
        None.

        """   
        self.nii_file = f1 
        nii1 = nib.load(f1)
        nii2 = nib.load(f2)

        im1in = nii1.get_fdata(dtype=np.float64)
        im2in = nii2.get_fdata(dtype=np.float64)
        self.snr_nema_2image(im1in, im2in)
        
    def snr_nema_multivol(self, nii_file):
        '''
        snr_nema_multivol
        alternative version of snr_nema, but for data where two images are
        from a single multi volume scan.

        Parameters
        ----------
        nii_file : STRING
            path to acquisition used in SNR calculation.

        Returns
        -------
        None.

        '''
        self.nii_file = nii_file
        nii = nib.load(nii_file)
        img = nii.get_fdata(dtype=np.float64)
        
        im1in = img[:,:,:,0]
        im2in = img[:,:,:,1]
        self.snr_nema_2image(im1in, im2in)
        
    def snr_nema_2image(self, im1in, im2in):
        '''
        snr_nema_2image(im1in,
                        im2in)
        Calculate SNR using the NEMA 2 image method: Two repeat images and 
        estimate SNR from the mean and difference images.

        Parameters
        ----------
        im1in : numpy memmap, generated by nibabel get_fdata
            first image in SNR calc.
        im2in : numpy memmap, generated by nibabel get_fdata
            second image in SNR calc.

        Returns
        -------
        None.

        '''
        if im1in.shape[0] != im2in.shape[0]:
            return None
        if im1in.shape[1] != im2in.shape[1]:
            return None
        if im1in.ndim != im2in.ndim:
            return None
        self.im_shape = im1in.shape        
        if self.im_shape[-1] != 1: #is a 3D nii
            # take midpoint in 
            mid_slice = int(np.floor( self.im_shape[-1] / 2 ))          
        else:
            # it's 2D
            mid_slice = 0
        self.im1 = im1in[:,:,mid_slice]
        self.im2 = im2in[:,:,mid_slice]
        
        self.locate_phantom(self.im1)

        # NEMA recommends using mask that's 75% of the area of the visible 
        # phantom, so that's sqrt(0.75) of the radius
        mask_radius = np.power(0.75,0.5) * self.ph_radius
        phantom_mask = self.mask_phantom(mask_radius)
 
    
        # NEMA Method 1: Use two images, estimate noise from subtraction
        im1_im2 = 0.5 * (self.im1 + self.im2)
        mean_signal = np.mean(im1_im2[phantom_mask==True])
        # subtraction for noise calculation
        self.im_subtract = self.im1 - self.im2
        im_sub_mask = self.im_subtract[phantom_mask==True] #unravelled
        std_signal = np.std(im_sub_mask)
        self.snr = np.power(2, 0.5) * mean_signal / std_signal
        self.mean_signal = mean_signal
        print('NEMA SNR (method 1: 2 image subtraction)       {:.6} '.format(self.snr))
        
        # NEMA Method 4 variant: Use single image, estimate noise from outside phantom
        #   this includes most of the region outside the phantom, so will include
        #   ghosting.
        noise_radius = 1.2 * self.ph_radius # estimated 20% larger than phantom
        noise_mask = self.mask_empty(noise_radius)       
        im_noise = self.im2[noise_mask==True] # use image 2 - less flow in phantom.
        sd_noise = np.std(im_noise)
        self.snr_bgd_ghost = 0.66 * np.mean(self.im1[phantom_mask==True]) / sd_noise
        print('     SNR (ref background including ghosting)   {:.6} '.format(self.snr_bgd_ghost))

        # NEMA Method 4 
        mask_corners = self.mask_corners(corner_length=20)
        im_corners = self.im2[mask_corners==True]
        sd_corners = np.std(im_corners)
        self.snr_bgd = 0.66 * np.mean(self.im1[phantom_mask==True]) / sd_corners
        print('NEMA SNR (ref corners)                         {:.6} '.format(self.snr_bgd))
       
        # Report the ghosting as the stdev of the whole region outside the phantom
        # (which should contain ghosting)
        # relative to the ghosting in the corners (which should be free of ghosting)
        self.ghost = sd_noise / sd_corners
        
        # plot       
        fig = plt.figure(figsize=(10,5))
        fig.suptitle(self.nii_file)
        fig.text(0,0.85, str("SNR_NEMA   {:.2f}".format(self.snr)))
        fig.text(0,0.8, str("SNR_bgd       {:.2f}".format(self.snr_bgd)))
        fig.text(0,0.75, str("Ghosting       {:.2f}".format(self.ghost)))

        ax = fig.subplots(2,2)
        ax[0][0].imshow(self.im1,cmap='jet')
        ax[0][1].imshow(self.im2,cmap='jet')
        ax[1][0].imshow(self.im_subtract ,cmap='jet')       
        ax[1][1].plot(self.im1[self.im_centre[1],:])  # dim0 in blue
        ax[1][1].plot(self.im1[:,self.im_centre[0]])  # dim1 in orange
        #ax[1][1].plot(im_sub_mask)
        
        # write report data, if required
        if self.write_report:
            self.output_root = self.nii_file.split('/')[-1].split('.nii')[0]
            self.output_png = self.output_root + '.png'
            fig.savefig(os.path.join(self.report_path, self.output_png))
            
        txt_fname = os.path.join(self.report_path, 'quickqc_'+self.output_root+'.dat')
        with open(txt_fname, 'w') as ff:
            ff.write('File,SNR_NEMA,SNR_background,Ghosting\n')
            ff.write( self.output_root \
                     +",{0:.2f}".format(self.snr)\
                     +",{0:.2f}".format(self.snr_bgd)\
                     +",{0:.2f}".format(self.ghost) )

        
    def uniformity_nema(self, f1):
        """        
        Principles:
            get data, 
            use same helper functions as snr to mask image.  Without masking, 
               there is a bimodal distribution of near zeros (outside phantom)
               and signal in the range ~1500 (inside phantom)
            descriptive statistics on the uniformity.
               Here a median probably makes more sense a head coil produces a 
               highly asymmetric distribution.  Median = value which has half
               of other values below it and half of other values above it.  
               It's the 50% value in pd.describe
            some plots (2D? thresholded?)
            

        Parameters
        ----------
        f1 : string
            Path to nii file for uniformity analysis.

        Returns
        -------
        None.

        """
        nii = nib.load(f1)
        im_in = nii.get_fdata(dtype=np.float64)        
        if im_in.shape[-1] != 1: #is a 3D nii
            # take midpoint in 
            mid_slice = int(np.floor( im_in.shape[-1] / 2 ))          
        else:
            # it's 2D
            mid_slice = 0
        im_u = im_in[:,:,mid_slice]
        
        self.locate_phantom(im_u)
        mask_radius = np.power(0.75,0.5) * self.ph_radius
        mask = self.mask_phantom(mask_radius)
        sig_mean = np.mean(im_u[mask==True])
        uniform_pd = pd.DataFrame(columns=['pixel values'] \
                                  , data=im_u[mask==True])
        uniform_pd.hist(bins=100)     
        uniform_desc = uniform_pd.describe()
        sig_median = uniform_desc.loc['50%', 'pixel values'] 
        uniform_desc['pixel scaled'] = (uniform_desc['pixel values'] / sig_median) - 1
        print(uniform_desc)


# MultiVolQc is generic for checking multi-volume MR data (e.g. fmri, qmt)
class MultiVolQc:
    '''
    MultiVolQc class: data and methods for dealing with multi-volume MRI data
    '''
    def __init__(self,filename, in_vivo=False, run_report=False, report_path=False):
        self.in_nii_file = os.path.basename(filename)
        self.in_nii_file_root = self.in_nii_file.split('.')[0]
        self.nii_path = os.path.dirname(filename)
        self.is_fmri = False
        self.is_diffusion = False
        self.in_vivo = in_vivo
        self.sfnr = 0
        # load the data here
        self.nii_load()
        self.basic_stats()
        if report_path == False:
            self.report_path = os.path.join(self.nii_path, self.in_nii_file_root +'_report')
        else:
            self.report_path = report_path
        if not os.path.exists(self.report_path):
            os.mkdir(self.report_path)
        if run_report:
            self.create_report()

    def nii_load(self):
        '''
        MultiVolQc.nii_load()

        Load nifti file using nibabel load.
        Populate nii_img (nibabel) and vol_data (numpy array)
        
        '''        
        self.nii_img = nib.load(os.path.join(self.nii_path, self.in_nii_file))
        # this will load as a proxy image
        # need to be careful with the datatypes better to load as float16 as this 
        # matches the dicom data and reduces memory usage, but some numpy functions
        # require float64
        self.vol_data = self.nii_img.get_fdata(dtype=np.float16)
        # reverse order of axes, so that we have [vol, slice, phase, read]
        self.vol_data = np.transpose(self.vol_data)
        self.shape = self.vol_data.shape
        self.affine = self.nii_img.affine
        if self.vol_data.ndim > 3:
            self.is_multi_volume = True
            self.n_vols = self.shape[0]
        else:
            self.is_multi_volume = False
            self.n_vols = 1
            
    def basic_stats(self, savenii=False, lthresh=0.1):
        '''
        MultiVolQc.basic_stats( 
            savenii=False
            )
        
        Calculate mean signal, stdev for each voxel and a simple mask
        '''
        # calculate volume mean, stdev and sfnr
        self.vol_mean = np.mean(self.vol_data,0, dtype=np.float64)
        self.vol_stdev = np.std(self.vol_data,0, dtype=np.float64)
        # 25% mask too high for 7T.  Reduce to 10%
        self.vol_mask = threshold_vol(self.vol_mean, True, lthresh)  
        if savenii:
            # create nifti, using same affine transform as original
            nii_mean = nib.Nifti1Image(self.vol_mean, self.affine)
            nii_stdev = nib.Nifti1Image(self.vol_stdev, self.affine)
            nii_mask = nib.Nifti1Image(self.vol_mask, self.affine)
            nib.save(nii_mean, os.path.join(self.report_path, 'fmriqc_mean.nii'))
            nib.save(nii_stdev, os.path.join(self.report_path, 'fmriqc_stdev.nii'))
            nib.save(nii_mask, os.path.join(self.report_path, 'fmriqc_mask.nii'))

    def voi(self, box_size):
        '''
        MultiVolQc.voi(
            box_size
        )
        
        Parameters:
        ----------
            box_size: 3 elemet tuple of dimensions (slices, rows, cols)

        
        Returns
        -------
            mask(nslice,nrows,ncols): 1 for pixels in mask, NaN outside
        '''
        # midpoints of slice, row, columns
        # for even numbers, midpoint will be start of second half of data (zero indexing)
        mid_points =  [int(np.floor(dim/2)) for dim in self.shape[1:]]
        start = []
        end = []
        for ii in range(0,3):
            start.append(mid_points[ii]-int(np.floor(box_size[ii])/2))
            end.append(start[ii]+box_size[ii])
        mask = np.tile(np.nan,self.shape[1:])
        mask[start[0]:end[0], start[1]:end[1], start[2]:end[2]]=1
        return mask   
    
    def timeseries(self, mask=None, plot=False, savepng=False):
        '''
        MultiVolQc.timeseries(
            mask=None
            plot=False
            savepng=False
            )
        
        Calculate mean signal across whole volume or masked volume for each 
        timepoint in 4D dataset.
        
        Parameters:
        ----------
            mask: supply a mask with 1's for included voxels, NaNs elsewhere
                         otherwise the volume_mask will be used
            plot:  generate a plot
            savepng: save plot
        
        Returns
        -------
        sig_timeseries

        '''
        
        if not np.any(mask):
            # mask at % of peak voxel intensity of mean image (set in basic_stats)
            mask = self.vol_mask 
        # keep volume data as float16 to hep=lp with memory problems
        masked_data = np.array(self.vol_data*mask, dtype=np.float16)
        # but polyfit needs float64 (but data smaller now)
        # and np.mean needs float64 too
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")        
            sig_timeseries = np.array( \
                np.nanmean(np.nanmean(np.nanmean(masked_data, axis=3, dtype=np.float64), axis=2, dtype=np.float64), axis=1, dtype=np.float64))
        if plot:
            fig = plt.figure(figsize=(30/2.56, 20/2.56))
            ax = fig.subplots(1,1)
            ax.plot(sig_timeseries)
            ax.set_xlabel('Volume No.')
            ax.set_ylabel('Mean signal')
            if savepng:
                fig.savefig(os.path.join(self.report_path, 'timeseries.png'))
        return sig_timeseries
                   
    def slice_time_plot(self, save_png=False):
         '''
         MultiVolQc.slice_time_plot(save_png=False)

         Plot mean signal from each slice for all time points (no masking)

         Parameters:
         save_fig:  Boolean.  Save figure as image
         
         '''
         slice_time = np.mean(np.mean(self.vol_data, axis=3), axis=2).T
         slice_mean = np.mean(np.mean(np.mean(self.vol_data, axis=3), axis=2),axis=0).T 
         slice_mean = np.tile(slice_mean,[self.shape[0],1]).T
         slice_time = slice_time - slice_mean

         fig = plt.figure(figsize=(30/2.5, 20/2.5))
         ax = fig.subplots()
         ax.set_title('Mean signal per slice, volume')
         ax.set_xlabel('Volume No.')
         ax.set_ylabel('Slice No.')
         im = ax.imshow(slice_time, cmap='plasma')
         if save_png:
             fig.savefig(os.path.join(self.report_path,'slice_time.png'))
                 

# methods specific to fmri
class FmriQc(MultiVolQc):
    '''
    methods specific to fmri (inherited from mriqc)
    '''
    def __init__(self,filename, in_vivo=False, run_report=True, report_path=False):
        MultiVolQc.__init__(self,filename, in_vivo, run_report, report_path)
        self.is_fmri = True
        
    def remove_dummies(self, dummies):
        '''
        In some conditions additional volumes need to be removed, above those
        removed by the sequence.  Yellow phantom on 7T is one case
        '''
        self.vol_data = self.vol_data[dummies:, :, :, :]

    def drift_correct(self, correct=False, mask=False, plot=False, savepng=False):
        '''
        FmriQc.drift_correct(
            correct=False
            plot=False
            savepng=False
            )
        Parameters
        ----------
        correct: apply correction (True/False).  If true, this updates fmriqc.vol_data
                 to the drift corrected value, otherwise just returns the drift and optionally
                 plots output
        Returns
        -------
        drift: the peak-to-peak variation in the fit (as % of mean)
        None.
        '''
        st = self.timeseries(mask=mask)
        vol_no = np.arange(0, len(st))
        p = np.polyfit(vol_no, st, 2)
        fitplot = p[0]*vol_no**2 + p[1]*vol_no + p[2]
        # correct linear and 2nd order term, but don't demean data (zeroth)
        corplot = np.array(p[0]*vol_no**2 + p[1]*vol_no, dtype=np.float16)
        resid=st-fitplot
        if plot:
            fig = plt.figure(figsize=(30/2.5, 20/2.5))
            ax1,ax2=fig.subplots(2,1)
            ax1.plot(vol_no, st, vol_no, fitplot)
            ax1.set_title('Pre drift correction, with fit')
            ax2.plot(vol_no, resid)
            ax2.set_title('Post drift correction')
            ax2.set_xlabel('Volume No.')

            if savepng:
                fig.savefig(os.path.join(self.report_path, 'drift_correct.png'))
                    
        drift = 100 * (np.amax(fitplot)-np.amin(fitplot)) /  np.mean(st)
        
        if correct:
            #cor = - np.tile(corplot[:,np.newaxis,np.newaxis,np.newaxis],
            #            [self.shape[1], self.shape[2], self.shape[3]] )
            #plt.plot(cor[:,18,32,32])
            #plt.set_title('cor')
            self.vol_data = self.vol_data - np.array(np.tile(corplot[:,np.newaxis,np.newaxis,np.newaxis], \
                        [self.shape[1], self.shape[2], self.shape[3]]), dtype=np.float16 )
            # if vol_data has been adjusted, need to recalculate mean, stdev
            self.basic_stats()   
    
        return(drift)
        
    def calc_sfnr(self, mask=None, plot=False, savepng=False, savenii=False):
        '''
        FmriQc.calc_sfnr(
            mask=False
            plot=True
            savepng=False
            savenii=False
            )
        calculate timeseries mean, standard deviation and signal to
                    fluctuation noise sfnr (Glover)
        Parameters:
            mask - 3D np array of booleans defining the mask, without a supplied
                   mask, the self.vol_mask will be used (from basic_stats())
            savepng - save ortho_view of SFNR, if required
        Returns:
            sfnr - SFNR from masked voi (or thresholded volume)        
        Outputs:
        self.vol_mean = mean signal across timeseries
        self.vol_stdev = standard deviation of signal across timepoints
        self.vol_sfnr = signal to fluctuation noise sfnr (Glover)
        '''
        #  if volume sfnr required (i.e. no mask specified), calculate based on 
        #  threshold of mean volume
        if not np.any(mask):
            # mask at percentage of peak voxel intensity of mean image (from basic_stats)
            mask = self.vol_mask
        # mask is 1 for pixels within mask, NaN for those outside
        # plot the unmasked version, but report stats on the masked version
        self.vol_sfnr_whole = np.divide(self.vol_mean, self.vol_stdev, \
                             out=np.zeros_like(self.vol_mean), \
                             where=self.vol_stdev!=0)
        self.vol_mean = self.vol_mean * mask
        self.vol_stdev = self.vol_stdev * mask

        # deal with inf.  
        self.vol_sfnr = np.divide(self.vol_mean, self.vol_stdev, \
                             out=np.zeros_like(self.vol_mean), \
                             where=self.vol_stdev!=0)
        if plot==True:
            ortho_view(self.vol_sfnr_whole, title='SFNR', save_png=savepng, save_dir=self.report_path)

        # discard zero values (as these are probably from the mask)
        sfnr = np.nanmean(np.ravel(self.vol_sfnr))
        vmean = np.nanmean(np.ravel(self.vol_mean))
        vstd = np.nanmean(np.ravel(self.vol_stdev))
        print("SFNR={0:.2f}   Mean={1:.2f}   stdev={2:.2f}\n".format(sfnr, vmean, vstd))
        if savenii:
            # create nifti, using same affine transform as original
            nii_sfnr = nib.Nifti1Image(np.array(self.vol_sfnr, dtype=np.float64), self.affine)
            nib.save(nii_sfnr, os.path.join(self.report_path, 'fmriqc_sfnr.nii'))      
        return (sfnr, vmean, vstd)
        
    def create_report(self):
        # build the elements needed, in case not run already
        #self.basic_stats()
        sfnr_vol,mean_vol,sd_vol = self.calc_sfnr(plot=True, savepng=True)

        # histogram of all image values (4D)
        plot_histogram(self.vol_data,save_png=True, save_path=self.report_path)
        self.slice_time_plot(True)
        
        if self.in_vivo:
            drift = self.drift_correct(correct=True,mask=False, plot=True, savepng=True)
#            t_series = self.timeseries(mask=None, plot=False, savepng=True)
        else:
            voi_mask = np.array(self.voi((10,20,20)), dtype=np.float16)
            self.remove_dummies(20)
            drift = self.drift_correct(correct=True,mask=voi_mask, plot=True, savepng=True)
#            t_series = self.timeseries(mask=voi_mask, plot=False, savepng=True)

        if not self.in_vivo:
            # running calc_sfnr with mask overwrites values in self
            sfnr_voi,mean_voi,sd_voi = self.calc_sfnr(voi_mask, plot=False, savenii=False)           
        else:
            sfnr_voi = 0
            mean_voi = 0
            sd_voi = 0  

        base_fname = self.in_nii_file.split('.')[0] # base filename from nii file
        # deal with XA long year format
        if base_fname[0:2]=='20':
            base_fname=base_fname[2:]
        html_fname = os.path.join(self.report_path, 'fmriqc_'+base_fname+'.html')
        with open(html_fname, 'w') as f:
            f.write('<!doctype=html><title>fMRI QC</title>\n<body>\n<p>fMRI QC Report</p>\n')     
            f.write('<table><tr><td>Image Dimensions</td><td>' + str(self.shape) + "</td></tr>\n")
            f.write('<tr><td>SFNR (volume) </td><td>' + "{0:.2f}".format(sfnr_vol) + "</td></tr>\n")
            f.write('<tr><td>Mean (volume) </td><td>' + "{0:.2f}".format(mean_vol) + "</td></tr>\n")
            f.write('<tr><td>stdev(volume) </td><td>' + "{0:.2f}".format(sd_vol) + "</td></tr>\n")
            f.write('<tr><td>SFNR (VOI) </td><td>' + "{0:.2f}".format(sfnr_voi) + "</td></tr>\n")
            f.write('<tr><td>Mean (VOI) </td><td>' + "{0:.2f}".format(mean_voi) + "</td></tr>\n")
            f.write('<tr><td>stdev(VOI) </td><td>' + "{0:.2f}".format(sd_voi) + "</td></tr>\n")
            f.write('<tr><td>Drift (%)</td><td>' + "{0:.2f}".format(drift) + "</td></tr>\n")
            f.write('</table>\n')
            pp = os.path.join('SFNR.png')
            f.write('<img src="' + pp + '"><br>\n')
            pp = os.path.join('pixel_histogram.png')
            f.write('<img src="' + pp + '"><br>\n')
            pp = os.path.join('drift_correct.png')      
            f.write('<img src="' + pp + '"><br>\n')
            pp = os.path.join('slice_time.png')      
            f.write('<img src="' + pp + '"><br>\n')
            f.write('</body>\n')
            
        txt_fname = os.path.join(self.report_path, 'fmriqc_'+base_fname+'.dat')
        with open(txt_fname, 'w') as ff:
            ff.write('File,sfnr_vol,mean_vol,sd_vol,sfnr_voi,mean_voi,sd_voi,drift\n')
            ff.write( base_fname \
                     +",{0:.2f}".format(sfnr_vol)\
                     +",{0:.2f}".format(mean_vol)\
                     +",{0:.2f}".format(sd_vol)\
                     +",{0:.2f}".format(sfnr_voi)\
                     +",{0:.2f}".format(mean_voi)\
                     +",{0:.2f}".format(sd_voi)\
                     +",{0:.2f}\n".format(drift)  )
        print("sfnr_volume =  {0:.2f}\n".format(sfnr_vol)\
            + "mean_volume =  {0:.2f}\n".format(mean_vol)\
            + "std_volume  =  {0:.2f}\n".format(sd_vol)\
            + "sfnr_VOI    =  {0:.2f}\n".format(sfnr_voi)\
            + "mean_VOI    =  {0:.2f}\n".format(mean_voi)\
            + "std_VOI     =  {0:.2f}\n".format(sd_voi)\
            + "drift (%)   =  {0:.2f}".format(drift) )

class MultiVolDiffusion(MultiVolQc):
    '''
    MultiVolDiffusion:
        Minimalistic diffusion class as a template for DiffusionQc

    '''
    def __init__(self,filename, in_vivo=True, run_report=False):
        MultiVolQc.__init__(self,filename, in_vivo=True, run_report=False)
        self.is_diffusion = True
        self.n_directions = self.read_bval()
        self.read_bvec()

    def read_bval(self):
        with open(os.path.join(self.nii_path,self.in_nii_file_root+'.bval')) as f:
            bvals_str=f.readline().replace('\n','').split(' ')
            bval_fl = []
            [ bval_fl.append(float(s)) for s in bvals_str ]
            self.bval = np.array(bval_fl)
            self.is_b0 = self.bval==True
            return self.bval.size
            
    def read_bvec(self):
        '''
        returns a bvec np array in the format bvec[axis][volume]
        '''        
        with open(os.path.join(self.nii_path,self.in_nii_file_root+'.bvec')) as f:
            bvec_fl = []
            for i in range(0,3):
                bvec_str=f.readline().replace('\n','').split(' ')
                [ bvec_fl.append(float(s)) for s in bvec_str ]
            self.bvec = np.array(bvec_fl).reshape([3,self.n_directions])   
    
class DiffusionQc(MultiVolDiffusion):
    '''
    DiffusionQc:
        DiffusionQc class for dealing with CUBRIC diffusion qc protocol
    
    '''
    def __init__(self,filename, in_vivo=True, run_report=False, nax=None, \
                 namp=None, nrep=None, nstart=None, nend=None, b0_in_loop=False):
        '''
        Parameters
        ----------
        filename : 
            nii_file 
        
        in_vivo:
            is in vivo
            
        run_report:
            generate html report
        
        nax:
            number of gradient axes in DiffusionQc acquisition (X,Y,Z)
        
        namp:
            number of gradient amplitudes in loop (may or may not include
            b0, depending on b0_in_loop)
            
        nrep:
            number of repeats, to help with SNR
            
        b0_in_loop: 
            True/False.  Are b0s in the main loop or at start/end?
            
        n_start: 
            vols at start (if outside the main loop)
            
        n_end:  
            vols at end (if outside the main loop)
            n_amp:



        Returns
        -------
        None.
     
        '''
        
        MultiVolDiffusion.__init__(self,filename, in_vivo=True, run_report=False)
        self.n_axes = nax
        self.n_amp = namp
        self.n_rep = nrep
        self.b0_in_loop = b0_in_loop
        self.n_start = nstart
        self.n_end = nend
    
    def prep_axis_amp_rep(self):
        '''
        prep_axis_amp_rep(
            )

        Prepare vol_data for gradient uniformity analysis. Use this function
        for vol_data which is in the format: Axis (slowest varying), amplitude, 
        repeats (fastest varying) i.e. 
        b0 (n_start)
        x, amp1 (n_rep)
        x, amp2 (n_rep)
        ...
        y, amp1 (n_rep)
        ...
        b0 (n_end)

        Returns
        -------
        None.

        '''
        
        b0 = self.vol_data[0:self.n_start,:,:,:]
        print(b0.shape)
        # number of vols in a block of gradient axis acquisitions is n_rep*n_amp
        # start geneating a set of views into original data
        x_start = self.n_start
        y_start = self.n_start + self.n_amp*self.n_rep
        z_start = y_start + self.n_amp*self.n_rep
        x_vol = self.vol_data[x_start:y_start]
        y_vol = self.vol_data[y_start:z_start]
        z_vol = self.vol_data[z_start:-self.n_end]

        self.x = DwGradAxis(b0,x_vol, self.vol_mask, self.voi([11,11,11]), self.n_rep, self.n_amp)
        self.y = DwGradAxis(b0,y_vol, self.vol_mask, self.voi([11,11,11]), self.n_rep, self.n_amp)
        self.x = DwGradAxis(b0,y_vol, self.vol_mask, self.voi([11,11,11]), self.n_rep, self.n_amp)
        
        plt.plot(x_vol[:,25,32,32])
        plt.title('X')
    
    def prep_rep_amp_axis(self):
        '''
        prep_rep_amp_axis(
            )

        Prepare vol_data for gradient uniformity analysis. Use this function
        for vol_data which is in the format: Repeats (slowest varying, biggest loop), 
        amplitude, axis (fastest varying) 
        i.e. 
        x, amp1, rep1
        y, amp1, rep1, 
        z, amp1, rep1
        x, amp2, rep1
        y, amp2, rep1
        ...
        x, amp1, rep2
        ...
        x, ampN, repM
        y, ampN, repM
        z, ampN, repM

        Returns
        -------
        None.

        '''
    
class DwGradAxis:
    '''
    Formatted data and methods for inspecting DiffusionQc data for a single gradient axis

    '''
    def __init__(self, b0_vols, grad_axis_vols, mask, voi_mask, n_rep, n_amp):
        '''
 
        Parameters
        ----------
        b0s_vols : numpy nd array
            An array of 3D vols, each of which is the non-dw image 

        grad_axis_vols : numpy nd array
            An array of 3D vols, each of which is the dw image corresponding to the dw signal 
            for gradient amplitude given by grad_amp.  Needs n_amp 3D vols
    
        grad_amp : list 
            Gradient amplitudes corresponding to n_amp measurements
            
        n_amp : int
            number of amplitude measurements

        Returns
        -------
        None.

        '''
        # check to see whether there are reps in b0
        # assign to b0 member of DeGradAxis
        if b0_vols.ndim > 3:
            self.b0 = np.mean(b0_vols,0)
        
        self.mask = mask
        self.voi_mask = voi_mask
        self.n_amp = n_amp

        if n_rep > 1:
            g1 = grad_axis_vols.reshape([n_amp,n_rep, \
                                          grad_axis_vols.shape[-3], \
                                          grad_axis_vols.shape[-2], \
                                          grad_axis_vols.shape[-1]])
            # mean across the n_reps dimension 
            self.diff_weighted = np.mean(g1,1)
        

    def g_uniformity(self):
        '''
        g_uniformity(
               )

        Params:
                diff_weighted    :  diffusion weighted (np array)
                non_diff_weighted:  b0 (np array)
        Returns:
                sqrt(ln(s0/sdiff)), i.e. the signal ratio linearised in G
                
        Some of the requirements of the protocol:
                The centre of the image has to align with scanner isocentre, as this is taken as the
                reference location for 'true' gradient.  An 11x11x11 voxel, centred on isocentre
                is taken as the reference position.

                Equation being solved here is 
                   --------------        
                -\/  ln(s0/sdiff)  = k(1+e) G_ideal
                where
                G_ideal is the ideal gradient (i.e. that requested)
                k is a constant containing timing terms, effective diffusion, gamma
                e is the required error term (as a fraction of G_ideal)
                
                k is estimated by assuming the gradient at isocentre is ideal

        '''
        print('s0_diff diff_weighted')
        print([self.b0.shape,self.diff_weighted.shape])
        
        s0_sdiff = np.multiply(np.divide(self.b0, self.diff_weighted), \
                        self.mask)
        sq_ln_sig = np.sqrt(np.log(s0_sdiff))

        # kG is the ideal (k * G) term - approximated by sq_ln_sig as isocentre
        kG_voi = np.multiply(sq_ln_sig, self.voi_mask)
        g_iso = np.nanmean(kG_voi)
        g_uniformity_vol = np.divide(sq_ln_sig,g_iso)
        print(g_uniformity_vol.shape)
        g_stdev = np.nanstd(g_uniformity_vol)
        for i in range(self.n_amp):
            ortho_view(g_uniformity_vol[i])
        plot_histogram(g_uniformity_vol)
        # need to deal with the fact that multiple volumes are reported here
        return(g_iso, g_stdev, g_uniformity_vol)
    
class SpikeQc(MultiVolQc):
    '''
    Methods for analysis of spike noise check data (single slice, multi vol EPI)
    call by creating a SpikeQc instance (passing filename), then calling
    spike_check()
    
    '''
    
    def spike_check(self, plot_title='spike'):
        '''
        SpikeQc.spike_check (
            
        )
        
        Params:
            
        Returns:
        '''
        
        slice_mean = self.timeseries(mask=None, plot=False,savepng=True)
        # indexes so slice, sorted in ascending slice_mean
        # require this to be filtered so that it's abouve 3std, then used
        # as fancy index to vol_data
        # sort ascending as default - want descending
        slice_mean_sort_idx = np.argsort(slice_mean)[::-1]
        # sorted, descending
        slice_mean_sorted = slice_mean[slice_mean_sort_idx]
        
        # spike slices - those with a mean signal > 3 stdev above the slice mean
        all_slice_mean = np.mean(slice_mean)
        all_slice_std = np.std(slice_mean)
        # sorted
        is_spike = slice_mean_sorted > all_slice_mean + 3*all_slice_std
        # only want the >3SD ones
        spike_idx = slice_mean_sort_idx[is_spike]
        # spike_slices contains image data with spikes, spikiest first
        self.spike_slices = self.vol_data[spike_idx,:,:,:]

        # mean sig plot, with outliers
        fig1,ax1 = plt.subplots(1,2)
        v=range(len(slice_mean))
        ax1[0].plot(v,slice_mean)
        ax1[0].scatter(spike_idx,slice_mean[spike_idx],marker='o', color='red')
        ax1[0].set_xlabel('Volume No.')
        ax1[0].set_ylabel('Mean signal')
        ax1[0].set_title(plot_title)

        #plot histogram
        plot_histogram(slice_mean, ax=ax1[1])   
        fig1.savefig(os.path.join(self.report_path,plot_title+'_spike_stats'))

        # plot slice images
        plot_row_max=5
        plot_col_max=5
        fig = plt.figure(figsize=(30/2.5, 30/2.5))
        ax = fig.subplots(plot_row_max,plot_col_max)
        for row in range(plot_row_max):
            for col in range(plot_col_max):
                ax[row][col].set_axis_off()
        if spike_idx.size > 15:
            plot_sl=15
        else:
            plot_sl=spike_idx.size
        used_rows = np.floor_divide(plot_sl,plot_col_max)
        for sl in range(plot_sl):
            rr = np.floor_divide(sl,plot_col_max)
            cc = np.mod(sl,plot_col_max)
            ax[rr][cc].imshow(self.spike_slices[sl,0,:,:], cmap='hsv')
        fig.savefig(os.path.join(self.report_path,plot_title+'_spike_images'))

            
class QcOverview():
    '''
    class QcOverview:
        Generate summary of QC data over a time period
        Requires an path_to_reports with output 
        *.dat files summarising outputs for a qc run.
        Uses os.walk(), so pointing to e.g. 
            /cubric/collab/108_QA/QA7T/fmriqc_glover/
        will navigate to the proc directory to find report *.dat files
        
    '''
    def __init__(self,system,path_to_reports, email_summary=False):
        '''
        Parameters
        ----------
        path_to_reports : string
            Path to top level where FmriQc reports are saved

        Requires the path_to_reports directory to have a 'summary' subdirectory for analysis.

        Returns
        -------
        None.
        '''
        self.dat_files = []
        self.system = system
        self.path_to_reports = path_to_reports 
        print('\nChecking for new ' + system + ' fmriqc analysis to add to summary')
        for root,dirs,files in os.walk(path_to_reports):
            for ff in files:
                if '.dat' in ff:
                    self.dat_files.append(os.path.join(root,ff))  
        run_analysis = self.check_analysis_required()
       
        # this launches the summary analysis 
        if run_analysis:
            print('Running summary analysis on' + path_to_reports)
            self.dat_to_pandas()
            self.plots()
            if email_summary:
                self.email_results()
        
        def dat_to_pandas(self):
            '''
            this is void for base class
            '''
            return False
        
        def plots(self):
            '''
            this is void for base class
            '''
            return False            
            
    def check_analysis_required(self):
        '''
        check_analysis_required:
            check whether new QA data has arrived since last analysis
            work out the time of;
            - last qa run
            - last time summary was run on qa data
            
        Returns
        -------
        new_qa = True/False. True if there is new QA data since last analysis

        '''
        
        # find latest download date
        t_data=self.find_latest_file(os.path.join(self.path_to_reports, 'proc' ) )
        print('latest QA dataset downloaded ' + str(time.ctime(t_data)))
        t_analysis = self.find_latest_file(os.path.join(self.path_to_reports, 'summary'))
        print('latest QA dataset analysed    ' + str(time.ctime(t_analysis)))
        if t_data > t_analysis:
            print('New QA data since last analysis')
            new_qa = True
        else:
            print('No new QA data')
            new_qa = False
        return new_qa
            
            
        # find lastst summary date
        
    def find_latest_file(self, check_path):
        '''
        navigate check_path to get lastest file time (as in since start 
        of time epoch)
        '''
        
        file_list = os.listdir(check_path)
        t_latest = 0
        for f in file_list:
            t = os.path.getctime(os.path.join(check_path, f))
            if t > t_latest:
                t_latest = t
        return t_latest
    
    def email_results(self):
        email_address = 'd6e057f1.cf.onmicrosoft.com@emea.teams.ms' # 3TW
        
        # Issue #2 - this hangs on wl026
        subprocess.run(['mailx','-s','fmriqc '+self.system, \
                        '-a', os.path.join(self.path_to_reports, 'summary','fmriqc_latest.png'), \
                        '-a', os.path.join(self.path_to_reports, 'summary','fmriqc_latest.txt'), \
                        email_address])
    
class FmriQcOverview(QcOverview):
    '''
    class FmriQcOverview.  Inherits generic methods for finding .dat files
    from QcOverview, but has specific methods for plotting and summarising 
    data
    '''
    def dat_to_pandas(self):
        '''
        Read in dat files (from self.dat_files) into pandas dataframe

        Returns
        -------
        None.

        '''
        data_str = []
        for df in self.dat_files:
            with open(df) as f:
                title=f.readline().replace('\n','').split(',')
                data_str.append(f.readline().replace('\n','').split(','))

        self.oview_qc=pd.DataFrame(data_str)
        self.oview_qc.columns = title
        self.oview_qc['sfnr_vol'] = self.oview_qc['sfnr_vol'].astype('float')
        self.oview_qc['mean_vol'] = self.oview_qc['mean_vol'].astype('float')
        self.oview_qc['sd_vol'] = self.oview_qc['sd_vol'].astype('float')
        self.oview_qc['sfnr_voi'] = self.oview_qc['sfnr_voi'].astype('float')
        self.oview_qc['mean_voi'] = self.oview_qc['mean_voi'].astype('float')
        self.oview_qc['sd_voi'] = self.oview_qc['sd_voi'].astype('float')
        self.oview_qc['drift'] = self.oview_qc['drift'].astype('float')
#        self.oview_qc.index=self.oview_qc['File']
#        self.oview_qc.drop('File', axis=1)
        date_str = self.oview_qc['File'].str.split('-', expand=True)[0]
        date_pd = date_str.str.split('_', expand=True)
        # deal with annoying XA patientid update.  Add 2000 if year is 2 digit.
        # year_temp could be YY (VB, VD, VE) or YYYY (XA)
        date_pd = date_pd.rename(columns={0:'year_temp',1:'month',2:'day'})
        date_pd = date_pd.astype('int16')
        def four_digit_year(yr):
            if yr < 100:
                return(yr + 2000)
            else:
                return(yr)
        date_pd['year']=date_pd['year_temp'].apply(four_digit_year)
        date_pd.pop('year_temp')
        self.oview_qc['date']=pd.to_datetime(date_pd)
        self.oview_qc['scanner'] = self.oview_qc['File'].str.split('-', expand=True)[3].str.split('_', expand=True)[8]
        self.oview_qc=self.oview_qc.rename(columns={'sfnr_vol':'sfnr_volume', 'mean_vol':'mean_volume', 'sd_vol':'sd_volume'})    
        self.oview_qc.tail()
        self.oview_qc_short = self.oview_qc
        self.oview_qc_short.index = self.oview_qc_short['date']
        self.oview_qc_short=self.oview_qc_short.drop(columns=['date','File'])
        self.oview_qc_short=self.oview_qc_short.sort_index(ascending=False)
        self.oview_qc_short=self.oview_qc_short.iloc[0:10]
        self.oview_qc_short.to_csv(os.path.join(self.path_to_reports, 'summary' \
                        , 'fmriqc_'+str(dt.datetime.now().isoformat()[:-7]).replace(':','-')+'.txt'))
        self.oview_qc_short.to_csv(os.path.join(self.path_to_reports, 'summary' \
                        , self.system+'_fmriqc_latest.txt'))
        
        
    def plots(self):
        fig = plt.figure(figsize=(18,6))
        axes = fig.subplots(1,3)
        self.oview_qc.plot(x='date', y=['sfnr_volume', 'sfnr_voi']\
                           , ax=axes[0], title = self.system+' EPI SFNR', grid=True, marker='o')   
        self.oview_qc.plot(x='date', y=['mean_volume', 'mean_voi']\
                           , ax=axes[1], title=self.system+' EPI Mean Signal', grid=True, marker='o')
        self.oview_qc.plot(x='date', y=['drift']\
                           , ax=axes[2], title=self.system+' EPI Drift (%)', grid=True, marker='o')
        fig.savefig(os.path.join(self.path_to_reports, 'summary'\
                        , 'fmriqc_'+str(dt.datetime.now().isoformat()[:-7]).replace(':','-')))
        fig.savefig(os.path.join(self.path_to_reports, 'summary'\
                        , self.system+'_fmriqc_latest.png'))
  
            
class SpikeQcOverview(QcOverview):
    '''
    Placeholder for spike analysis
    '''            
            

def threshold_vol(vol, by_fraction, threshold):
    '''
    mriqc.threshold_vol(
       by_fraction
       threshold
    )

    Params:
    by_fraction:  bool.
       True = threshold_val is a fraction of max value in image
       False = threshold_val is absolute pixel intensity
    threshold: float/int
       np.nan for all values below threshold

    Returns:
        np array with mask of values above threshold
    '''

    max_pixel = np.amax(np.amax(np.amax(vol)))
    # mask has to be a copy of vol, to prevent the original volume being overwritten by the mask
    mask = np.array(vol)
    pixel_threshold = max_pixel * threshold
    mask[mask<pixel_threshold] = np.nan
    mask[mask>=pixel_threshold] = 1
    return mask

def ortho_view(vol, title='image', save_png=False, save_dir='.', im_scale=[]):
    '''
    mriqc.ortho_view(
        vol, 
        title='image', 
        savepng=False,
        save_dir='.'
        im_scale=[vmin,vmax]
        )
    Params:
    vol: 3D numpy array to display.  Can't be multi_volume
    title: plot title
    save_png: save output as png
    save_dir: save to directory
    im_scale: set image scale limits to [vmin, vmax]
    Show middle slice in three orthogonal views
    '''

    if vol.ndim != 3:
        print('vol has ' + str(vol.ndim) + ' volumes')
        print('Must be a 3D numpy array')
        return

    mid_slice = [int(np.floor(dim_len/2)) for dim_len in vol.shape]
    orth = []
    vmax = -10000
    vmin = 10000
    orth.append(vol[:,:,mid_slice[2]])
    orth.append(vol[:,mid_slice[1],:])
    orth.append(vol[mid_slice[0],:,:])

    # get max, min values from orth slices if not supplied by im_scale
    if len(im_scale)==2:
        vmin = im_scale[0]
        vmax = im_scale[1]
    else:
        for sl in orth:
            if np.amax(sl[~np.isnan(sl)]) > vmax:
                vmax = np.amax(sl[~np.isnan(sl)])
            if np.amin(sl[~np.isnan(sl)]) < vmin:
                vmin = np.amin(sl[~np.isnan(sl)])
    
    fig = plt.figure(figsize=(30/2.56, 20/2.56))
    ax = fig.subplots(2,2)
    im1=ax[0,0].imshow(orth[0], origin='lower', cmap='gray', vmax=vmax, vmin=vmin)
    im2=ax[0,1].imshow(orth[1], origin='lower', cmap='gray', vmax=vmax, vmin=vmin)
    im3=ax[1,0].imshow(orth[2], origin='lower', cmap='gray', vmax=vmax, vmin=vmin)
    for r in [0,1]:
        for c in [0,1]:
             ax[r,c].axis('off')
    ax[1,1].set_title(title)
    fig.colorbar(im3) #place the colormap in bottom left
    if save_png:
        fig.savefig(os.path.join(save_dir, title + '.png'))
 
    return mid_slice

def plot_histogram(vol, save_png=False, save_path = '.', ax = None):
    '''
    plot_histogram(vol, save_png=False, save_path='.'):

    Parameters:
    vol = 3D or 4D np array
    save_png = save output as png
    save_path = directory to save png

    
    Plots a histogram of all pixel values, across all values in 3D or 4D
    numpy array.  Zero values discarded 
    Saves histogram in a png file
    '''

    hist = np.histogram(vol[vol>0],100)
    bins = hist[1][1:]  # bins includes both ends, take the higher value
    vals = hist[0]
    fig = plt.figure(figsize=(30/2.56, 20/2.56))
    if ax == None:
        ax = fig.subplots()
    ax.plot(bins, vals)
    ax.set_xlabel('Pixel Value')
    ax.set_ylabel('Number of pixels')
    ax.set_title('Image histogram')
    if save_png:
        fig.savefig(os.path.join(save_path, 'pixel_histogram.png'))
    return ax
        
def prep_data(scanner,root='/cubric/collab/108_QA'):
    '''
    mriqc.prep_data(
        scanner: one of { '3TE', '3TW', '7T', 'connectom'}
        root='/cubric/collab/108_QA': the path_root of QC data
        )
        
    Prepares QC data for analysis.  Requires the following folder structure:
      .
      ├── scanner
      │   ├── fmriqc_2023
      │   ├── fmriqc_glover
      │   │   ├── proc     
      │   │   │   ├── session1
      │   │   │   │   ├──  nifti1
      │   │   │   │   └──  report1          
      │   │   │   └── session2
      │   │   │       ├──  nifti2
      │   │   │       └──  report2  
      │   │   └── summary                  
      │   ├── nifti
      │   │   ├── session1
      │   │   │   ├──  nifti1
      │   │   │   └──  nifti2          
      │   │   └── session2 
      │   │       ├──  nifti1
      │   │       └──  nifti2                   
      │   ├── raw
      │   │   ├──  session1
      │   │   └──  session2    
      │   ├── spike_2017
      │   └── spike_2023
  
    Session directories should be downloaded by Download Images on XNAT with 
    'simplify directory struture' OFF.
    e.g.
    ├── 23_05_24-15_12_01-DST-1_3_12_2_1107_5_2_34_18984
│   └── scans
│       ├── 1-localizer_gradiso
│       ├── 2-gre_2dms_coilcheck
│       ├── 3-gre_2dms_coilcheck
│       ├── 4-spike_readLR
│       ├── 5-spike_readAP
│       ├── 6-spike_readHF
│       ├── 7-QC_MB_GRE_EPI_FA15_Tx220V_SBRef
│       └── 8-QC_MB_GRE_EPI_FA15_Tx220V

    
    Does the following:
        1) look through all root/scanner/raw for a session folder
        2) check whether session has nifti data in root/scanner/nifti
        3) if yes - skip
           if no  - run dcm2niix (linux only)
        4) move niftis to root/scanner/qctype/proc/nifti or    
    
    
    Parameters
    ----------
    
        
    Returns
    -------
    None.

    '''
    import shutil
    import subprocess
    
    all_scanners = ['3TE', '3TW', '7T', 'connectom']
    if scanner not in all_scanners:
        print('!! ERROR: ' + scanner +' is not one of ')
        print([s for s in all_scanners])
        return
    
    dt = os.path.join(root,scanner) #top directory
    dr = os.path.join(dt, 'raw') #raw dir
    dn = os.path.join(dt, 'nifti') #nifti dir
    sessr = []         #list of sess
    for f in os.listdir(dr):
        if '-1_3_12_2_1107_5_2_' in f:
            sessr.append(f)
    
    sessn = []
    for ff in os.listdir(dn):
        if '-1_3_12_2_1107_5_2_' in ff:
            sessn.append(ff)
    
    # check all session directories
    for s in sessr:
        print('\n\n\n      ###  ' + s + '  ###\n\n\n')
        # check that there isn't a matching session in nifti directory
        # if not, do dicom conversion, and move resulting files to nifti dir
        if s not in sessn:
            print('No niftis for ' + s + '. Running dcm2niix.')
            spath = os.path.join(dr,s)
            npath = os.path.join(dn,s)
            subprocess.run(['echo','/cubric/software/bin/dcm2niix', '-f', '%i_%s_%d', spath])
            subprocess.run(['/cubric/software/bin/dcm2niix', '-f', '%i_%s_%d', spath])
            os.mkdir(npath)
            fn = [ f for f in os.listdir(spath) if 'nii' in f ]
            print(fn)
            for niftifile in fn:
                shutil.move(os.path.join(spath,niftifile), npath)       
            fj = [ f2 for f2 in os.listdir(spath) if 'json' in f2 ]
            for jsonfile in fj:
                shutil.move(os.path.join(spath, jsonfile), npath)
        else:
            print('Niftis already exist for '+ s + '. Skipping.')
           
    
    # need to update nifti list, as may have been added to.
    sessniinew = []
    for ff in os.listdir(dn):
        if '-1_3_12_2_1107_5_2_' in ff:
            sessniinew.append(ff)
            
    # then check nifti directory and move to relevant qc directory
    p_fmriqc_2023 = os.path.join(root,scanner,'fmriqc_2023/proc')
    p_fmriqc_glover= os.path.join(root, scanner, 'fmriqc_glover/proc')
    p_spike_2023 = os.path.join(root, scanner, 'spike_2023/proc')
    p_spike_2017 = os.path.join(root, scanner, 'spike_2017/proc')
    
    for sn in sessniinew:
        print('Sorting niftis for ' + sn)
        snpath = os.path.join(dn,sn)
        niiscans = os.listdir(snpath)
        for nii in niiscans:
            # fmriqc_2023
            if ('QC_MB_GRE_EPI' in nii) and ('.nii' in nii) and ('SBRef' not in nii):
                print(nii + ' >> ' + p_fmriqc_2023)
                shutil.move(os.path.join(snpath,nii), p_fmriqc_2023)
            # spike_2023
            if ('spike_read' in nii) and ('.nii' in nii):
                print(nii + ' >> ' + p_spike_2023)
                shutil.move(os.path.join(snpath,nii), p_spike_2023)
            if ('GloverGSQAP' in nii) and ('.nii' in nii):
                print(nii + ' >> ' + p_fmriqc_glover) 
                shutil.move(os.path.join(snpath,nii), p_fmriqc_glover)
            if ('EPIspike' in nii) and ('.nii' in nii):
                print(nii + ' >> ' + p_spike_2017)
                shutil.move(os.path.join(snpath,nii), p_spike_2017)

    



    
