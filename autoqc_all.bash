#!/bin/bash
# Wrapper script for cron to run all qc setup commands in one file, and deal with logfiles 

log_file='/cubric/collab/108_QA/logs/autoqc.log'

date > ${log_file}
echo "----------------   autoqc_status.py (pre)  ----------------" >> ${log_file}
/home/sapje1/code/python_mrobjects/qc/autoqc_status.py >> ${log_file}

echo "----------------   autoqc_xnat_fetch.py    ----------------" >> ${log_file}
/home/sapje1/code/python_mrobjects/qc/autoqc_xnat_fetch.py >> ${log_file}

echo "----------------   autoqc_process.py       ----------------" >> ${log_file}
/home/sapje1/code/python_mrobjects/qc/autoqc_process.py >> ${log_file}

echo "----------------   autoqc_status.py (post) ----------------" >> ${log_file}
/home/sapje1/code/python_mrobjects/qc/autoqc_status.py >> ${log_file}

/home/sapje1/code/python_mrobjects/qc/autoqc_email.bash

mailx -s "autoqc.log" evansj31@cardiff.ac.uk  < ${log_file}

