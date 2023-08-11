# -*- coding: utf-8 -*-
"""
caliber_distortion analysis
CJE 10/8/2023
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
    
    def flip_hyperfine_directions(self):
        '''
        flip the directions reported by the Caliber output to match the physical directions of the Hyperfine;
        translation is:
            Caliber S->I : Hyperfine A->P
            Caliber L->R : Hyperfine L->R (or R->L?)
            Caliber A->P : Hyperfine I->S
        '''
        self.distortion.rename( \
                      columns={'L/R Expected': 'Expected LR' \
                      , 'A/P Expected': 'Expected SI' \
                      , 'S/I Expected': 'Expected AP' \
                      , 'Measured L/R': 'Measured LR' \
                      , 'Measured A/P': 'Measured SI' \
                      , 'Measured S/I': 'Measured AP' \
                      , 'L/R Distortion': 'Distortion LR' \
                      , 'A/P Distortion': 'Distortion SI' \
                      , 'S/I Distortion': 'Distortion AP' \
                      }, inplace=True )  
    
    def plot_distortion(self):
        fig = plt.figure(figsize=(18,12))
        fig.suptitle(self.study_title)
        axes = fig.subplots(3,3)
        d=self.distortion
        d.plot(x='Expected LR', y='Distortion LR',ylabel='Distortion LR (mm)',style='xb',ax=axes[0][0], legend=False)
        d.plot(x='Expected LR', y='Distortion AP',ylabel='Distortion AP (mm)',style='xg',ax=axes[0][1], legend=False)
        d.plot(x='Expected LR', y='Distortion SI',ylabel='Distortion SI (mm)',style='xg',ax=axes[0][2], legend=False)
    
        d.plot(x='Expected AP', y='Distortion LR',ylabel='Distortion LR (mm)',style='vg',ax=axes[1][0], legend=False)
        d.plot(x='Expected AP', y='Distortion AP',ylabel='Distortion AP (mm)',style='vb',ax=axes[1][1], legend=False)
        d.plot(x='Expected AP', y='Distortion SI',ylabel='Distortion SI (mm)',style='vg',ax=axes[1][2], legend=False)
        
        d.plot(x='Expected SI', y='Distortion LR',ylabel='Distortion LR (mm)',style='og',ax=axes[2][0], legend=False)
        d.plot(x='Expected SI', y='Distortion AP',ylabel='Distortion AP (mm)',style='og',ax=axes[2][1], legend=False)
        d.plot(x='Expected SI', y='Distortion SI',ylabel='Distortion SI (mm)',style='ob',ax=axes[2][2], legend=False)
    
    def plot_measured(self):
        fig = plt.figure(figsize=(18,12))
        fig.suptitle(self.study_title)
        axes = fig.subplots(1,3) 
        d=self.distortion

        d.plot(x='Expected LR', y='Measured LR',ylabel='Measured LR (mm)',style='xb', \
               markersize=10, ax=axes[0], legend=False)
        d.plot(x='Expected AP', y='Measured AP',ylabel='Measured AP (mm)',style='1g', markersize=10,ax=axes[1], legend=False)
        d.plot(x='Expected SI', y='Measured SI',ylabel='Measured SI (mm)',style='.r', markersize=10,ax=axes[2], legend=False)
        
        x_ax = np.array([-50,50])
    
        [si, ap, lr] = self.grad_calibration()
    
        # fit  = x * linear_term  + constant
        si_fit_plot = x_ax * si[0]   +   si[1]
        ap_fit_plot = x_ax * ap[0]   +   ap[1]
        lr_fit_plot = x_ax * lr[0]   +   lr[1]
    
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
    
    
    def grad_calibration(self):
        si_fit = np.polyfit(self.distortion.loc[:,'Expected SI'], self.distortion.loc[:,'Measured SI'], deg=1)
        ap_fit = np.polyfit(self.distortion.loc[:,'Expected AP'], self.distortion.loc[:,'Measured AP'], deg=1)
        lr_fit = np.polyfit(self.distortion.loc[:,'Expected LR'], self.distortion.loc[:,'Measured LR'], deg=1)
        return [si_fit, ap_fit, lr_fit]
    
    def plot_stats(self):
        fig = plt.figure(figsize=(18,6))
        fig.suptitle(self.study_title)
        axis = fig.subplots(1,1) 
        self.distortion.loc[:,['Dist From Expected','Distortion LR', 'Distortion AP', 'Distortion SI']].plot.box(ylabel='Distortion (mm)', ax=axis, title='Summary Stats')    
    
    def analyse_hyperfine(self):
        self.read_caliber_distortion()
        self.sort_caliber_distortion()
        self.flip_hyperfine_directions()
        self.plot_distortion()
        self.plot_measured()
        self.plot_stats()
    
        [si, ap, lr] = self.grad_calibration()
        si_grad_scaling = (si[0]-1) * 100
        lr_grad_scaling = (lr[0]-1) * 100
        ap_grad_scaling = (ap[0]-1) * 100
        grad_scaling = [lr_grad_scaling, ap_grad_scaling, si_grad_scaling]
        b0_iso_shift = [lr[1], ap[1], si[1]]
        
        print('Gradient scaling LR: {:.3f} %'.format(lr_grad_scaling))
        print('Gradient scaling AP: {:.3f} %'.format(ap_grad_scaling))
        print('Gradient scaling SI: {:.3f} %'.format(si_grad_scaling))
        print('Iso shift (B0)   LR: {:.3f} mm'.format(lr[1]))
        print('Iso shift (B0)   AP: {:.3f} mm'.format(ap[1]))
        print('Iso shift (B0)   SI: {:.3f} mm'.format(si[1]))


