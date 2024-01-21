''' 
 autoqc_fetch.py
 wrapper script to call the various qc functions and modules to automate
 the qc analysis
 CJE Sept 2023
'''

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import xnat_fetch_qc as xnqc

# get new qc data from xnat
xnqc.update_xnat_new()
xnqc.xnat_download()
xnqc.data_unzip()
