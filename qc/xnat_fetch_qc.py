#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xnat_fetchqc.py

Collect qc data from xnat project 108 using xnatpy
Created on Wed Aug 16 14:37:20 2023

@author: sapje1
"""

import xnat
import subprocess
import os

data_path = '/cubric/collab/108_QA'
download_list = os.path.join(data_path, 'xnat_download_list.txt')

qcsubj = { 'QA7T' : 'XNAT_S06014',
           'QA3TM' : 'XNAT_S06051',
           'QA3TE' : 'XNAT_S06097',
           'QA3TW' : 'XNAT_S06091'
        }

def update_download_list():
    #regenerate list of downloaded datasets
    exp_done = []
    print('Checking for downloaded data in ' + data_path)
    for qd in qcsubj.keys():
        qdir = os.path.join(data_path,qd,'raw')
        list1 = os.listdir(qdir)
    
        for l in list1:
            if 'XNAT_E' in l and '.zip' not in l:
                exp_done.append(l)
    print('Found ' + str(len(exp_done)) + ' downloaded datasets:')
    print(exp_done)
    with open(download_list,'w') as f:
        for line in exp_done:
            f.write(line + '\n')

def xnat_download():
    # phase 2: download new zip files from XNAT
    with open(download_list, 'r') as f:
        tmp = f.readlines()
    exp_downloaded = []
    for line in tmp:
        exp_downloaded.append(line.replace('\n',''))
        
    # if no authentication provided, xnat will look to .netrc file for authentication
    print('\nChecking for new data on XNAT')
    
    with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
        subj = session.projects['108_QA2023'].subjects
        for s in subj:
            for (q1, q2) in qcsubj.items():
                if q2 in s:
                    subj_id = subj[s].label
                    qc_exp=subj[s].experiments
                    print('Found ' + str(len(qc_exp)) + ' sessions on XNAT for ' + subj_id)
                    # e is all qc experiments for a given subject (=scanner)
                    exp_to_download = []
                    for e in qc_exp:
                        if e not in exp_downloaded:
                            raw_path = os.path.join(data_path, subj_id, 'raw', e + '.zip')
                            print('Session ' + e + ' is new.')
                            exp_to_download.append(e)
                    print(str(len(exp_to_download)) + ' new session will be downloaded')
                    for ed in exp_to_download:
                        print('Downloading ' + e + ' to ' + raw_path)
                        qc_exp[e].download(raw_path)                                               
    session.disconnect()

def data_unzip():
    # phase 3: unzip data to /cubric/collab/108_QA/SCANNER/raw/XNATID/SESSIONID
    print('\nChecking for XNAT zip files in ' + data_path)
    for q1, q2 in qcsubj.items():
        p = os.path.join(data_path,q1,'raw')
        fls = os.listdir(p)
        for ff in fls:
            if 'XNAT' in ff[0:4] and '.zip' in ff[-4:]:
                zip_dir = os.path.join(p,ff[:-4])
                ff_path = os.path.join(p,ff)
                print('Unzipping ' + ff_path + ' to ' + zip_dir)
                subprocess.run(['unzip', '-q', '-d', zip_dir, ff_path])

#update_download_list()
#xnat_download()
data_unzip()