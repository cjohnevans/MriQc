#!/bin/bash
# Wrapper script for cron to run all qc setup commands in one file, and deal with logfiles 

log_file='/cubric/collab/108_QA/logs/autoqc.log'

date > ${log_file}
echo "----------------   autoqc_status.py (pre)  ----------------" >> ${log_file}
/home/sapje1/code/MriQc/autoqc_status.py >> ${log_file}

echo "----------------   autoqc_process.py       ----------------" >> ${log_file}
/home/sapje1/code/MriQc/autoqc_process.py >> ${log_file}

echo "----------------   autoqc_status.py (post) ----------------" >> ${log_file}
/home/sapje1/code/MriQc/autoqc_status.py >> ${log_file}

/home/sapje1/code/MriQc/autoqc_email.bash

#mailx -s "autoqc.log" evansj31@cardiff.ac.uk  < ${log_file}
cat ${log_file} | mailx -r mriqc@cardiff.ac.uk -s "autoqc.log" 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
