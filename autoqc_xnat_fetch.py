#!/home/sapje1/miniconda2/envs/mri/bin/python3

''' 
 autoqc_xnat_fetch.py

 Fetch available data from XNAT, based on information in /cubric/collab/108_QA2023 files;
     - xnat_download_done.txt
     - xnat_new.txt

 !!! To refresh these, need to run autoqc_status.py !!!

 ORDER:
    autoqc_status.py
    autoqc_xnat_fetch.py
    autoqc_process.py


 CJE Jan 2024
'''

import sys
sys.path.append('/home/sapje1/code/MriQc')
import xnat_fetch_qc as xnqc

#  currently not using this - require manual check of autoqc_status prior to auto_qc_fetch
# update file with list of local downloads
#xnqc.update_downloaded()
# get new qc data from xnat
#xnqc.update_xnat_new()

xnqc.xnat_download()

