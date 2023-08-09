# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

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