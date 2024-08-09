#!/home/sapje1/miniconda2/envs/mri/bin/python3

''' 
 autoqc_process.py

  Does the following: 
  - Unzip data from XNAT 
  - Convert unzipped data to nifti 
  - Process QC data
  - Summarise data


 ORDER:
    autoqc_status.py
    autoqc_xnat_fetch.py
    autoqc_process.py
 
 CJE Jan 2024
'''

import sys
sys.path.append('/home/sapje1/code/MriQc')
import xnat_fetch_qc as xnqc
import mriqc

# get new qc data from xnat
xnqc.data_unzip(unzip=True, remove_invalid_file=True)
xnqc.nifti_convert()

# process QC data
xnqc.proc_qc()


