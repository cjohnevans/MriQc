import nibabel as nib
import numpy as np
import os
import matplotlib.pyplot as plt

# MultiVolQc is generic for checking multi-volume MR data (e.g. fmri, qmt)
class MultiVolQc:
    '''
    MultiVolQc class: data and methods for dealing with multi-volume MRI data
    '''
    def __init__(self,path,name, in_vivo=True):
        self.in_nii_file = name
        self.nii_path = path
        self.in_vivo = False # default to phantom
        self.sfnr = 0
        # load the data here
        self.nii_load()
        self.basic_stats()
        self.report_path = os.path.join(self.nii_path, 'report')
        if not os.path.exists(self.report_path):
            os.mkdir(self.report_path)

    def nii_load(self):
        '''
        MultiVolQc.nii_load()

        Load nifti file using nibabel load.
        Populate nii_img (nibabel) and vol_data (numpy array)
        
        '''        
        self.nii_img = nib.load(os.path.join(self.nii_path, self.in_nii_file))
        # this will load as a proxy image
        self.vol_data = self.nii_img.get_fdata()
        # reverse order of axes, so that we have [vol, slice, phase, read]
        self.vol_data = np.transpose(self.vol_data)
        self.shape = self.vol_data.shape
        self.affine = self.nii_img.affine
        if self.vol_data.ndim > 3:
            self.is_multi_volume = True
            self.n_vols = self.shape[-1]
        else:
            self.is_multi_volume = False
            self.n_vols = 1
            
    def basic_stats(self, savenii=False):
        '''
        MultiVolQc.basic_stats( 
            savenii=False
            )
        
        Calculate mean signal, stdev for each voxel and a simple mask
        '''
        # calculate volume mean, stdev and sfnr
        self.vol_mean = np.mean(self.vol_data,0)
        self.vol_stdev = np.std(self.vol_data,0)
        self.vol_mask = threshold_vol(self.vol_mean, True, 0.25)  
        if savenii:
            # create nifti, using same affine transform as original
            nii_mean = nib.Nifti1Image(self.vol_mean, self.affine)
            nii_stdev = nib.Nifti1Image(self.vol_stdev, self.affine)
            nii_mask = nib.Nifti1Image(self.vol_mask, self.affine)
            nib.save(nii_mean, os.path.join(self.nii_path, 'fmriqc_mean.nii'))
            nib.save(nii_stdev, os.path.join(self.nii_path, 'fmriqc_stdev.nii'))
            nib.save(nii_mask, os.path.join(self.nii_path, 'fmriqc_mask.nii'))
    
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
            # mask at 25% of peak voxel intensity of mean image
            mask = self.vol_mask 
        
        masked_data = self.vol_data*mask
        sig_timeseries = np.nanmean(np.nanmean(np.nanmean(masked_data, axis=3), axis=2), axis=1)
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
    def drift_correct(self, correct=False, plot=False):
        '''
        FmriQc.drift_correct(
            correct=False
            plot=False
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
        st = self.timeseries()
        vol_no = np.arange(0, len(st))
        p = np.polyfit(vol_no, st, 2)
        fitplot = p[0]*vol_no**2 + p[1]*vol_no + p[2]
        resid=st-fitplot
        if plot:
            fig=plt.figure()
            ax1,ax2=fig.subplots(2,1)
            ax1.plot(vol_no, st, vol_no, fitplot)
            ax2.plot(vol_no, resid)

        drift = 100 * (np.amax(fitplot)-np.amin(fitplot)) /  np.mean(st)
        
        if correct:
            self.vol_data = self.vol_data - np.tile(fitplot[:,np.newaxis,np.newaxis,np.newaxis],
                        [self.shape[1], self.shape[2], self.shape[3]] )
            
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
        self.basic_stats()
        #  if volume sfnr required (i.e. no mask specified), calculate based on 
        #  threshold of mean volume
        if not np.any(mask):
            # mask at 25% of peak voxel intensity of mean image
            mask = self.vol_mask
        # mask is 1 for pixels within mask, NaN for those outside
        self.vol_mean = self.vol_mean * mask
        self.vol_stdev = self.vol_stdev * mask

        # deal with inf.  
        self.vol_sfnr = np.divide(self.vol_mean, self.vol_stdev, \
                             out=np.zeros_like(self.vol_mean), \
                             where=self.vol_stdev!=0)
        if plot==True:
            ortho_view(self.vol_sfnr, title='SFNR', save_png=savepng, save_dir=self.report_path)

        # discard zero values (as these are probably from the mask)
        sfnr = np.nanmean(np.ravel(self.vol_sfnr))
        print(sfnr)        
        if savenii:
            # create nifti, using same affine transform as original
            nii_sfnr = nib.Nifti1Image(self.vol_sfnr, self.affine)
            nib.save(nii_sfnr, os.path.join(self.nii_path, 'fmriqc_sfnr.nii'))      
        return sfnr
        
    def create_report(self):
        # build the elements needed, in case not run already
        self.basic_stats() 
        sfnr_vol=self.calc_sfnr(plot=True, savepng=True)

        # histogram of all image values (4D)
        plot_histogram(self.vol_data,save_png=True, save_path=self.report_path)
        self.slice_time_plot(True)
        
        if self.in_vivo:
            t_series = self.timeseries(mask=None, plot=False)
        else:
            voi_mask = self.voi((10,20,20))
            t_series = self.timeseries(mask=voi_mask, plot=False)

        html_fname = os.path.join(self.report_path, 'fmriqc_report.html')
        with open(html_fname, 'w') as f:
            f.write('<!doctype=html><title>fMRI QC</title>\n<body>\n<p>fMRI QC Report</p>\n')     
            f.write('<table><tr><td>Image Dimensions</td><td>' + str(self.shape) + "</td></tr>\n")
            f.write('<tr><td>SFNR (volume) </td><td>' + "{0:.2f}".format(sfnr_vol) + "</td></tr>\n")
            if not self.in_vivo:
                sfnr_voi=self.calc_sfnr(voi_mask, plot=False)           
                f.write('<tr><td>SFNR (VOI) </td><td>' + "{0:.2f}".format(sfnr_voi) + "</td></tr>\n")
            f.write('<tr><td>Drift</td><td>' + "{0:.2f}".format(1.2345) + "</td></tr>\n")
            f.write('</table>\n')
            pp = os.path.join('SFNR.png')
            f.write('<img src="' + pp + '"><br>\n')
            pp = os.path.join('pixel_histogram.png')
            f.write('<img src="' + pp + '"><br>\n')
            pp = os.path.join('slice_time.png')      
            f.write('<img src="' + pp + '"><br>\n')
            f.write('</body>\n')

    def voi(self, box_size):
        '''
        FmriQc.voi(
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

def ortho_view(vol, title='image', save_png=False, save_dir='.'):
    '''
    mriqc.ortho_view(
        vol, 
        title='image', 
        savepng=False, 
        )
    Params:
    vol: 3D numpy array to display.  Can't be multi_volume
    title: plot title
    save_png: save output as png
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

    # get max, min values from orth slices
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

def plot_histogram(vol, save_png=False, save_path = '.'):
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
    ax = fig.subplots()
    ax.plot(bins, vals)
    ax.set_xlabel('Pixel Value')
    ax.set_ylabel('Number of pixels')
    ax.set_title('Image histogram')
    if save_png:
        fig.savefig(os.path.join(save_path, 'pixel_histogram.png'))


    
