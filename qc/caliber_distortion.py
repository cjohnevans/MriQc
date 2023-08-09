# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_caliber_distortion(dist_file):
    with open(dist_file) as f:
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
    return block

def sort_caliber_distortion(block):
    for b in block:
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
    return [study_param, series_param, distortion]

def flip_hyperfine_directions(distortion):
    '''
    flip the directions reported by the Caliber output to match the physical directions of the Hyperfine;
    translation is:
        Caliber S->I : Hyperfine A->P
        Caliber L->R : Hyperfine L->R (or R->L?)
        Caliber A->P : Hyperfine I->S
    '''
    dist2 = distortion.rename( \
                  columns={'L/R Expected': 'Expected LR' \
                  , 'A/P Expected': 'Expected SI' \
                  , 'S/I Expected': 'Expected AP' \
                  , 'Measured L/R': 'Measured LR' \
                  , 'Measured A/P': 'Measured SI' \
                  , 'Measured S/I': 'Measured AP' \
                  , 'L/R Distortion': 'Distortion LR' \
                  , 'A/P Distortion': 'Distortion SI' \
                  , 'S/I Distortion': 'Distortion AP' \
                          } )  
    return dist2

def plot_distortion(d, title):
    fig = plt.figure(figsize=(18,12))
    fig.suptitle(title)
    axes = fig.subplots(3,3) 
    d.plot(x='Expected LR', y='Distortion LR',ylabel='Distortion LR (mm)',style='xb',ax=axes[0][0], legend=False)
    d.plot(x='Expected LR', y='Distortion AP',ylabel='Distortion AP (mm)',style='xg',ax=axes[0][1], legend=False)
    d.plot(x='Expected LR', y='Distortion SI',ylabel='Distortion SI (mm)',style='xg',ax=axes[0][2], legend=False)

    d.plot(x='Expected AP', y='Distortion LR',ylabel='Distortion LR (mm)',style='vg',ax=axes[1][0], legend=False)
    d.plot(x='Expected AP', y='Distortion AP',ylabel='Distortion AP (mm)',style='vb',ax=axes[1][1], legend=False)
    d.plot(x='Expected AP', y='Distortion SI',ylabel='Distortion SI (mm)',style='vg',ax=axes[1][2], legend=False)
    
    d.plot(x='Expected SI', y='Distortion LR',ylabel='Distortion LR (mm)',style='og',ax=axes[2][0], legend=False)
    d.plot(x='Expected SI', y='Distortion AP',ylabel='Distortion AP (mm)',style='og',ax=axes[2][1], legend=False)
    d.plot(x='Expected SI', y='Distortion SI',ylabel='Distortion SI (mm)',style='ob',ax=axes[2][2], legend=False)

def plot_measured(d, title):
    fig = plt.figure(figsize=(18,12))
    fig.suptitle(title)
    axes = fig.subplots(1,3) 
    d.plot(x='Expected LR', y='Measured LR',ylabel='Measured LR (mm)',style='xb', \
           markersize=10, ax=axes[0], legend=False)
    d.plot(x='Expected AP', y='Measured AP',ylabel='Measured AP (mm)',style='1g', markersize=10,ax=axes[1], legend=False)
    d.plot(x='Expected SI', y='Measured SI',ylabel='Measured SI (mm)',style='.r', markersize=10,ax=axes[2], legend=False)
    
    x_ax = np.array([-50,50])

    [si, ap, lr] = grad_calibration(d)

    # fit  = x * linear_term  + constant
    si_fit = x_ax * si[0]   +   si[1]
    ap_fit = x_ax * ap[0]   +   ap[1]
    lr_fit = x_ax * lr[0]   +   lr[1]

#   plot nominal gradient and actual fit (first order only)
    axes[0].plot(x_ax, [-50, 50], 'k--')
    axes[0].plot(x_ax, lr_fit, 'b-')
    axes[1].plot(x_ax, [-50, 50], 'k--')
    axes[1].plot(x_ax, ap_fit, 'g-')
    axes[2].plot(x_ax, [-50, 50], 'k--')
    axes[2].plot(x_ax, si_fit, 'r-')
    for a in axes:
        t=np.arange(-50,51,10)
        a.set_yticks(t)
        a.set_xticks(t)
        a.grid(axis='both')


def grad_calibration(distortion):
    si_fit = np.polyfit(distortion.loc[:,'Expected SI'], distortion.loc[:,'Measured SI'], deg=1)
    ap_fit = np.polyfit(distortion.loc[:,'Expected AP'], distortion.loc[:,'Measured AP'], deg=1)
    lr_fit = np.polyfit(distortion.loc[:,'Expected LR'], distortion.loc[:,'Measured LR'], deg=1)
    return [si_fit, ap_fit, lr_fit]

def plot_stats(d, fig_title):
    fig = plt.figure(figsize=(18,6))
    fig.suptitle(fig_title)
    axis = fig.subplots(1,1) 
    d.loc[:,['Dist From Expected','Distortion LR', 'Distortion AP', 'Distortion SI']].plot.box(ylabel='Distortion (mm)', ax=axis, title='Summary Stats')
    

def analyse_hyperfine(dist_file):
    block = read_caliber_distortion(dist_file)
    [study, series, distortion] = sort_caliber_distortion(block)
    d2 = flip_hyperfine_directions(distortion)
    fig_title = study.loc['Study Description',0] + ': ' + series.loc['Sequence Name',0]
    plot_distortion(d2,fig_title)
    plot_measured(d2,fig_title)
    plot_stats(d2, fig_title)

    [si, ap, lr] = grad_calibration(d2)
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
          
    return [study, series, d2, grad_scaling]