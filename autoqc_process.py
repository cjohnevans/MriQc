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
    autoqc_process.py
 
 CJE Jan 2024
'''

import sys
sys.path.append('/home/sapje1/code/MriQc')
import mriqc_xnat as xnqc
import mriqc

qc_types = ['fmriqc', 'spikehead', 'spikebody']

for qc in qc_types:
    xnqc.update_downloaded(qc_type=qc)
    xnqc.update_xnat_new(qc_type=qc)
    xnqc.xnat_download(qc_type=qc)
    xnqc.data_unzip(qc_type=qc, unzip=True, remove_invalid_file=True)
    xnqc.nifti_convert(qc_type=qc)
    xnqc.proc_qc(qc_type=qc)

mriqc.FmriQcOverview('3TW','/cubric/collab/108_QA/QA3TW/fmriqc/', email_summary=False)
mriqc.FmriQcOverview('7T','/cubric/collab/108_QA/QA7T/fmriqc/', email_summary=False)
mriqc.FmriQcOverview('3TE','/cubric/collab/108_QA/QA3TE/fmriqc/', email_summary=False)
mriqc.FmriQcOverview('3TM','/cubric/collab/108_QA/QA3TM/fmriqc/', email_summary=False)
