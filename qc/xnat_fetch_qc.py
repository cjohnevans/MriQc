#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xnat_fetchqc.py

Collect qc data from xnat project 108 using xnatpy
Created on Wed Aug 16 14:37:20 2023

@author: sapje1
"""

import xnat
import sys
import os

data_path = '/cubric/collab/108_QA'
download_list = os.path.join(data_path, 'xnat_download_list.txt')
qcsubj = { 'QA7T' : 'XNAT_S06014',
           'QA3TM' : 'XNAT_S06051',
           'QA3TE' : 'XNAT_S06097',
           'QA3TW' : 'XNAT_S06091'
        }

with open(download_list) as f:
    tmp = f.readlines()
exp_downloaded = []
for line in tmp:
    exp_downloaded.append(line.replace('\n',''))

print(exp_downloaded)

# if no authentication provided, xnat will look to .netrc file for authentication
with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
    subj = session.projects['108_QA2023'].subjects


    for s in subj:
        for (q1, q2) in qcsubj.items():
            if q2 in s:
                subj_id = subj[s].label
                qc_exp=subj[s].experiments
                print('N sessions = ' + str(len(qc_exp)))
                for e in qc_exp:
                    if e in exp_downloaded:
                        print(e + ' downloaded')
                    else:
                        print(e)
                        raw_path = os.path.join(data_path, subj_id, 'raw', e + '.zip')
                        print(raw_path)
                        qc_exp[e].download(raw_path)
                        

session.disconnect()
