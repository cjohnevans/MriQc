import nibabel as nib
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import warnings

# MultiVolQc is generic for checking multi-volume MR data (e.g. fmri, qmt)
class MultiVolQc:
    '''
    MultiVolQc class: data and methods for dealing with multi-volume MRI data
    '''
    def __init__(self,filename, in_vivo=False, run_report=False):
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
        self.report_path = os.path.join(self.nii_path, self.in_nii_file_root +'_report')
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
    def __init__(self,filename, in_vivo=False, run_report=True):
        MultiVolQc.__init__(self,filename, in_vivo, run_report)
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
    
    '''
    
    def spike_check(self):
        '''
        SpikeQc.spike_check (
            
        )
        
        Params:
            
        Returns:
        '''
        
        slice_mean = self.timeseries(mask=None, plot=True,savepng=True)
        plot_histogram(slice_mean)
        
        # spike slices - those with a mean signal > 3 stdev above the slice mean
        all_slice_mean = np.mean(slice_mean)
        all_slice_std = np.std(slice_mean)
        spike_slices = slice_mean > all_slice_mean + 1.5*all_slice_std
        vv = np.arange(0,spike_slices.shape[0])
        spike_slice_idx = vv[spike_slices]
        self.spike_slices = self.vol_data[spike_slices,:,:,:]
        print(spike_slice_idx)
        plot_row_max=5
        plot_col_max=5
         
        fig = plt.figure(figsize=(30/2.5, 30/2.5))
        ax = fig.subplots(plot_row_max,plot_col_max)
        for row in range(plot_row_max):
            for col in range(plot_col_max):
                ax[row][col].set_axis_off()
        if spike_slice_idx.size > 15:
            plot_sl=15
        else:
            plot_sl=spike_slice_idx.size
        used_rows = np.floor_divide(plot_sl,plot_col_max)
        for sl in range(plot_sl):
            rr = np.floor_divide(sl,plot_col_max)
            cc = np.mod(sl,plot_col_max)
            ax[rr][cc].imshow(self.spike_slices[sl,0,:,:])
            
class FmriQcOverview():
    '''
    class FmriQcOverview:
        Generate summary of FmriQc data over a time period
        Requires an input path with output *_report directories containing
        summary *.dat files with sfnr etc outputs for a qc run.
    '''
    def __init__(self,path_to_reports):
        '''
        Parameters
        ----------
        path_to_reports : string
            Path to top level where FmriQc reports are saved 

        Returns
        -------
        None.
        '''
        self.dat_files = []
        for root,dirs,files in os.walk(path_to_reports):
            for ff in files:
                if '.dat' in ff:
                    self.dat_files.append(os.path.join(root,ff))                    
        self.dat_to_pandas()
        self.plots()
        
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
        date_pd = date_pd.rename(columns={0:'year_short',1:'month',2:'day'})
        date_pd = date_pd.astype('int16')
        date_pd['year'] = date_pd['year_short'].apply(lambda x: x+2000)
        date_pd.pop('year_short')
        self.oview_qc['date']=pd.to_datetime(date_pd)
        self.oview_qc['scanner'] = self.oview_qc['File'].str.split('-', expand=True)[3].str.split('_', expand=True)[8]
        self.oview_qc=self.oview_qc.rename(columns={'sfnr_vol':'sfnr_volume', 'mean_vol':'mean_volume', 'sd_vol':'sd_volume'})    
        self.oview_qc.tail()
        
    def plots(self):
        self.oview_qc.plot(x='date', y=['sfnr_volume', 'sfnr_voi'])   
        self.oview_qc.plot(x='date', y=['mean_volume', 'mean_voi'])
        self.oview_qc.plot(x='date', y=['sd_volume', 'sd_voi'])
        self.oview_qc.plot(x='date', y=['drift'])
          

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
        
        



    
