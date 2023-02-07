import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import os
import numpy as np
import matplotlib.pyplot as plt

if os.name == 'nt':
    nii_path = 'C:\\Users\\sapje1\\data\\RoutineQA_examples'
    nii_file = '20211111_113225GloverGSQAPs003a001.nii'
else:
    nii_path = '/home/sapje1/scratch_sapje1/fmriqc/251_alspac/'
    nii_file = 'alspacfmri.nii'

fmri = mriqc.fmriqc(nii_path,nii_file,in_vivo=False)

fmri.create_report()

t_series = fmri.timeseries(plot=True)
vol_no = np.arange(0,len(t_series))
t_series2 = t_series + 0.1*vol_no

p = np.polyfit(vol_no, t_series2, 2)
fitplot = p[0]*vol_no**2 + p[1]*vol_no + p[2]
resid = t_series2-fitplot

fig = plt.figure()
ax = fig.subplots()
ax.plot(vol_no, t_series2, vol_no, fitplot) 
fig2 = plt.figure()
ax2 = fig2.subplots()
ax2.plot(vol_no, resid)
