''' 
 autoqc_check.py
 
 check the status of the downloaded data
 
 CJE Sept 2023
'''

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import xnat_fetch_qc as xnqc

# get new qc data from xnat
xnqc.data_unzip(unzip=True, remove_invalid_file=True)
xnqc.nifti_convert()
