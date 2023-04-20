#!/bin/bash
source activate mri
/home/sapje1/miniconda2/envs/mri/bin/python3 /home/sapje1/code/python_mrobjects/qc/magmon7T_wl026.py | mailx -s "magmon7T" -a  /cubric/scanners/mri/7t/magnet_logs/cryostat_pressure.png -a /cubric/scanners/mri/7t/magnet_logs/cryostat_T.png 15f0963a.cf.onmicrosoft.com@uk.teams.ms 
