#!/home/sapje1/miniconda2/envs/mri/bin/python3

''' 
 autoqc_status.py
 
 check the status of;
     - the downloaded data
     - new data available on xnat
     - existing analyses
 
  Updates the download/available file lists
     - xnat_download_done.txt
     - xnat_new.txt


 ORDER:
    autoqc_status.py
    autoqc_xnat_fetch.py
    autoqc_process.py


 CJE Jan 2024
'''

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import xnat_fetch_qc as xnqc

# get list of downloaded datasets from /cubric/collab/108_QA2023
xnqc.update_downloaded()

# get new qc data from xnat
xnqc.update_xnat_new()

# check status of latest sessions
def read_last_qc(scanner, fname):
    with open(fname,'r') as ff:
        header=ff.readline()
        lastqc=ff.readline() 
    return lastqc

scannerlist = ['7T ', '3TM', '3TW','3TE']
filelist = ['/cubric/collab/108_QA/QA7T/fmriqc_glover/summary/7T_fmriqc_latest.txt', \
           '/cubric/collab/108_QA/QA3TM/fmriqc_glover/summary/3TM_fmriqc_latest.txt', \
           '/cubric/collab/108_QA/QA3TW/fmriqc_glover/summary/3TW_fmriqc_latest.txt', \
           '/cubric/collab/108_QA/QA3TE/fmriqc_glover/summary/3TE_fmriqc_latest.txt']

print('\nautoqc_status: Checking for latest QC analyses')
for ii in range(0,4):
    last_qc = read_last_qc(scannerlist[ii], filelist[ii])
    print(scannerlist[ii] + ' last QC: ' + last_qc[:-1])

