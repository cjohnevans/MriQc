''' 
 autoqc_proc.py
 wrapper script to call the various qc functions and modules to automate
 the qc analysis
 CJE Jan 2024
'''

import sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc
import xnat_fetch_qc as xnqc

xnqc.proc_qc()

# run fmriqc summary script
mriqc.FmriQcOverview('3TW','/cubric/collab/108_QA/QA3TW/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('7T','/cubric/collab/108_QA/QA7T/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('7T','/cubric/collab/108_QA/QA7T/fmriqc_MB/', email_summary=False)
mriqc.FmriQcOverview('3TE','/cubric/collab/108_QA/QA3TE/fmriqc_glover/', email_summary=False)
mriqc.FmriQcOverview('3TM','/cubric/collab/108_QA/QA3TM/fmriqc_glover/', email_summary=False)



