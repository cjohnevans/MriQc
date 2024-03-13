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
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import xnat_fetch_qc as xnqc
import mriqc

# get new qc data from xnat
xnqc.data_unzip(unzip=True, remove_invalid_file=True)
xnqc.nifti_convert()

# process QC data
xnqc.proc_qc()

# run fmriqc summary script
mriqc.FmriQcOverview('3TW','/cubric/collab/108_QA/QA3TW/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('7T','/cubric/collab/108_QA/QA7T/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('7T','/cubric/collab/108_QA/QA7T/fmriqc_MB/', email_summary=False)
mriqc.FmriQcOverview('3TE','/cubric/collab/108_QA/QA3TE/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('3TM','/cubric/collab/108_QA/QA3TM/fmriqc_glover/', email_summary=False)

mriqc.QuickQcOverview('3TW', '/cubric/collab/108_QA/QA3TW/quick_SNR_gre3D', email_summary=False)
mriqc.QuickQcOverview('7T', '/cubric/collab/108_QA/QA7T/quick_SNR_gre3D', email_summary=False)
mriqc.QuickQcOverview('3TE', '/cubric/collab/108_QA/QA3TE/quick_SNR_gre3D', email_summary=False)
mriqc.QuickQcOverview('3TM', '/cubric/collab/108_QA/QA3TM/quick_SNR_gre3D', email_summary=False)
