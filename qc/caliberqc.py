# -*- coding: utf-8 -*-
"""
caliber_distortion analysis
CJE 10/8/2023
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class GradientAxes():
    def __init__(self, lr, ap, si, name):
        self.lr=lr
        self.ap=ap
        self.si=si
        self.name = name
    def show(self):
        print(self.name)
        print('LR: {:.3f}'.format(self.lr))
        print('AP: {:.3f}'.format(self.ap))
        print('SI: {:.3f}'.format(self.si))


class CaliberDistortion():
    def __init__(self, dist_file):
        self.gradient_scaling = {}
        self.iso_shift = {}
        self.dist_file = dist_file

    def read_caliber_distortion(self):
        with open(self.dist_file) as f:
            templine=f.readline()
            # read data in blocks, for reading with pandas later
            block=[] # as list
            block_temp = []
            while templine: #until EOF
                if len(templine) > 1 :
                    tidytemp1=templine.split('\",\"')
                    tidytemp2 = [ t.replace('\"','').replace('\n','') for t in tidytemp1 ]
                    block_temp.append(tidytemp2)
                else: #end of block
                    if len(block_temp) > 0:    #ignore blank lines
                        block.append(block_temp) # block is a list of lists
                    block_temp = [] #reset
                templine=f.readline()
                
            # catch case where there's no blank line after a block at EOF
            if len(block_temp) > 0:
                block.append(block_temp) # block is a list of lists   
        self.block = block

    
    def sort_caliber_distortion(self):
        for b in self.block:
            title=b[0][0]
            if title=='Study Parameters':
                study_param=pd.DataFrame(b[1:], columns=b[0]).T
            if title=='Series Parameters':
                series_param=pd.DataFrame(b[1:], columns=b[0]).T
            if title=='VOI Statistics':
                distortion = pd.DataFrame(b[2:], columns=b[0])
                distortion['L/R Expected']=distortion['L/R Expected'].astype('float64')
                distortion['A/P Expected']=distortion['A/P Expected'].astype('float64')
                distortion['S/I Expected']=distortion['S/I Expected'].astype('float64')
                distortion['Measured L/R']=distortion['Measured L/R'].astype('float64')
                distortion['Measured A/P']=distortion['Measured A/P'].astype('float64')
                distortion['Measured S/I']=distortion['Measured S/I'].astype('float64')
                distortion['L/R Distortion']=distortion['L/R Distortion'].astype('float64')
                distortion['A/P Distortion']=distortion['A/P Distortion'].astype('float64')
                distortion['S/I Distortion']=distortion['S/I Distortion'].astype('float64')
                distortion['Center Dist vs Expected']=distortion['Center Dist vs Expected'].astype('float64')
                distortion['Dist From Expected']=distortion['Dist From Expected'].astype('float64')
        self.study_param =study_param
        self.series_param = series_param
        self.distortion = distortion
        self.study_title = study_param.loc['Study Description',0] + ': ' + series_param.loc['Sequence Name',0]
        self.Rmax = np.max(np.power(np.power(self.distortion['L/R Expected'],2) \
                             + np.power(self.distortion['A/P Expected'],2) \
                             + np.power(self.distortion['S/I Expected'],2), 0.5))
        self.dist_max = np.amax([np.amax(self.distortion['L/R Distortion']) \
                             , np.amax(self.distortion['A/P Distortion'] )  \
                             , np.amax(self.distortion['S/I Distortion']) ])
        self.dist_min = np.amin([np.amin(self.distortion['L/R Distortion']) \
                                 , np.amin(self.distortion['A/P Distortion'] )  \
                                 , np.amin(self.distortion['S/I Distortion']) ])
    
    def flip_hyperfine_directions(self):
        '''
        flip the directions reported by the Caliber output to match the physical directions of the Hyperfine;
        translation is:
            Caliber S->I : Hyperfine A->P
            Caliber L->R : Hyperfine L->R (or R->L?)
            Caliber A->P : Hyperfine I->S
        '''
        self.distortion.rename( \
                      columns={'L/R Expected': 'Expected LR (mm)' \
                      , 'A/P Expected': 'Expected SI (mm)' \
                      , 'S/I Expected': 'Expected AP (mm)' \
                      , 'Measured L/R': 'Measured LR (mm)' \
                      , 'Measured A/P': 'Measured SI (mm)' \
                      , 'Measured S/I': 'Measured AP (mm)' \
                      , 'L/R Distortion': 'Distortion LR (mm)' \
                      , 'A/P Distortion': 'Distortion SI (mm)' \
                      , 'S/I Distortion': 'Distortion AP (mm)' \
                      }, inplace=True )  
    
    def plot_distortion(self):
        fig = plt.figure(figsize=(18,12))
        fig.suptitle(self.study_title)
        axes = fig.subplots(3,3)
        d=self.distortion
        d.plot(x='Expected LR (mm)', y='Distortion LR (mm)',ylabel='Distortion LR (mm)',style='xb',ax=axes[0][0], legend=False)
        d.plot(x='Expected LR (mm)', y='Distortion AP (mm)',ylabel='Distortion AP (mm)',style='xg',ax=axes[0][1], legend=False)
        d.plot(x='Expected LR (mm)', y='Distortion SI (mm)',ylabel='Distortion SI (mm)',style='xg',ax=axes[0][2], legend=False)
    
        d.plot(x='Expected AP (mm)', y='Distortion LR (mm)',ylabel='Distortion LR (mm)',style='vg',ax=axes[1][0], legend=False)
        d.plot(x='Expected AP (mm)', y='Distortion AP (mm)',ylabel='Distortion AP (mm)',style='vb',ax=axes[1][1], legend=False)
        d.plot(x='Expected AP (mm)', y='Distortion SI (mm)',ylabel='Distortion SI (mm)',style='vg',ax=axes[1][2], legend=False)
        
        d.plot(x='Expected SI (mm)', y='Distortion LR (mm)',ylabel='Distortion LR (mm)',style='og',ax=axes[2][0], legend=False)
        d.plot(x='Expected SI (mm)', y='Distortion AP (mm)',ylabel='Distortion AP (mm)',style='og',ax=axes[2][1], legend=False)
        d.plot(x='Expected SI (mm)', y='Distortion SI (mm)',ylabel='Distortion SI (mm)',style='ob',ax=axes[2][2], legend=False)
        
        for a in range(len(axes)):
            for b in range(len(axes[a])): 
                t=np.arange(-10*np.ceil(self.Rmax/10),10*np.ceil(self.Rmax/10),10)
                axes[a][b].set_xticks(t)
                axes[a][b].grid(axis='both')           
    
    
    def plot_measured(self):
        fig = plt.figure(figsize=(18,12))
        fig.suptitle(self.study_title)
        axes = fig.subplots(1,3) 
        d=self.distortion

        d.plot(x='Expected LR (mm)', y='Measured LR (mm)',ylabel='Measured LR (mm)',style='xb', \
               markersize=10, ax=axes[0], legend=False)
        d.plot(x='Expected AP (mm)', y='Measured AP (mm)',ylabel='Measured AP (mm)',style='1g', markersize=10,ax=axes[1], legend=False)
        d.plot(x='Expected SI (mm)', y='Measured SI (mm)',ylabel='Measured SI (mm)',style='.r', markersize=10,ax=axes[2], legend=False)
        
        x_ax = np.array([-50,50])
    
        # fit  = x * linear_term  + constant
        si_fit_plot = x_ax * self.grad_scale.si  + self.grad_offset.si
        ap_fit_plot = x_ax * self.grad_scale.ap  + self.grad_offset.ap
        lr_fit_plot = x_ax * self.grad_scale.lr  + self.grad_offset.lr
    
    #   plot nominal gradient and actual fit (first order only)
        axes[0].plot(x_ax, [-50, 50], 'k--')
        axes[0].plot(x_ax, lr_fit_plot, 'b-')
        axes[1].plot(x_ax, [-50, 50], 'k--')
        axes[1].plot(x_ax, ap_fit_plot, 'g-')
        axes[2].plot(x_ax, [-50, 50], 'k--')
        axes[2].plot(x_ax, si_fit_plot, 'r-')
        for a in axes:
            t=np.arange(-50,51,10)
            a.set_yticks(t)
            a.set_xticks(t)
            a.grid(axis='both')
    
    
    def evaluate_gradient(self):
        si_fit = np.polyfit(self.distortion.loc[:,'Expected SI (mm)'], self.distortion.loc[:,'Measured SI (mm)'], deg=1)
        ap_fit = np.polyfit(self.distortion.loc[:,'Expected AP (mm)'], self.distortion.loc[:,'Measured AP (mm)'], deg=1)
        lr_fit = np.polyfit(self.distortion.loc[:,'Expected LR (mm)'], self.distortion.loc[:,'Measured LR (mm)'], deg=1)
        self.grad_scale = GradientAxes(lr_fit[0], ap_fit[0], si_fit[0], 'Gradient Linear Scaling Factor')
        self.grad_scale_pct=GradientAxes(  (lr_fit[0]-1) * 100 \
                            , (ap_fit[0]-1) * 100 \
                            , (si_fit[0]-1) * 100 \
                            , 'Gradient Linear Scale Error (%)')

        # iso_shift is the actual isocentre position (mm)
        # y=mx+c with y=0 ...
        self.grad_offset = GradientAxes ( lr_fit[1], ap_fit[1], si_fit[1], 'Grad fit constant offset')
        self.iso_shift = GradientAxes((- lr_fit[1] / lr_fit[0]) \
                                      ,(- ap_fit[1] / ap_fit[0]) \
                                      , (- si_fit[1] / si_fit[0]) \
                                      , 'Isocentre Shift (mm)')
        self.distortion['Distortion LR (%R)'] = 100 * self.distortion['Distortion LR (mm)'] / self.Rmax
        self.distortion['Distortion AP (%R)'] = 100 * self.distortion['Distortion AP (mm)'] / self.Rmax
        self.distortion['Distortion SI (%R)'] = 100 * self.distortion['Distortion SI (mm)'] / self.Rmax

        
        # here the nonlinearity is calculated as a difference from the predicted gradient, rather than
        # the fitted gradient, so it will include contributions from gradient linear scaling and non-orthogonality
        n_lr = self.distortion['Distortion LR (%R)'].describe()['max'] - \
                    self.distortion['Distortion LR (%R)'].describe()['min']
        n_ap = self.distortion['Distortion AP (%R)'].describe()['max'] - \
                    self.distortion['Distortion AP (%R)'].describe()['min']
        n_si = self.distortion['Distortion SI (%R)'].describe()['max'] - \
                    self.distortion['Distortion SI (%R)'].describe()['min']
        self.nonlinearity_peak = GradientAxes(n_lr, n_ap, n_si, 'Gradient Non-linearity (%, peak-to-peak)')
        r_lr = np.sqrt( np.mean ( np.power(self.distortion['Distortion LR (%R)'],2)))
        r_ap = np.sqrt( np.mean ( np.power(self.distortion['Distortion AP (%R)'],2)))     
        r_si = np.sqrt( np.mean ( np.power(self.distortion['Distortion SI (%R)'],2)))
        self.nonlinearity_rms = GradientAxes(r_lr, r_ap, r_si, 'Gradient Non-linearity (%, RMS)')

    def plot_stats(self):
        fig = plt.figure(figsize=(18,6))
        fig.suptitle(self.study_title)
        axis = fig.subplots(1,1) 
        self.distortion.loc[:,['Dist From Expected','Distortion LR (mm)', 'Distortion AP (mm)', 'Distortion SI (mm)']].plot.box(ylabel='Distortion (mm)', ax=axis, title='Summary Stats')    
    
    def analyse_hyperfine(self):
        self.read_caliber_distortion()
        self.sort_caliber_distortion()
        self.flip_hyperfine_directions()
        self.plot_distortion()
        self.evaluate_gradient()
        #self.plot_measured()
        self.plot_stats()
    
        self.grad_scale.show()
        self.grad_scale_pct.show()
        self.iso_shift.show()
        



