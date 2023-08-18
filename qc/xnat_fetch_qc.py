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
import os, shutil

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
        qdir = os.path.join(data_path,qd,'nifti')
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
            for (qc_subj, qc_subj_xn) in qcsubj.items():
                if qc_subj_xn in s:
                    subj_id = subj[s].label  # subj_id is QA3TM etc.. 
                    qc_exp = subj[s].experiments
                    print('Found ' + str(len(qc_exp)) + ' sessions on XNAT for ' + subj_id)
                    # e is all qc experiments for a given subject (=scanner)
                    exp_to_download = []
                    
                    # set up zip path for this scanner (subject)
                    dir_zip = os.path.join(data_path, subj_id, 'zip')
                    if not os.path.isdir(dir_zip):
                        os.mkdir(dir_zip)                     
                    
                    for e in qc_exp:
                        # check xnat experiement ids against those previously downloaded
                        if e not in exp_downloaded:
                            print('Session ' + e + ' is new.')
                            exp_to_download.append(e)
                            
                    print(str(len(exp_to_download)) + ' new session(s) will be downloaded')
                    for ed in exp_to_download:
                        zip_path = os.path.join(dir_zip, ed + '.zip')
                        print('Downloading ' + ed + ' to ' + zip_path)
                        qc_exp[ed].download(zip_path)                                               
    session.disconnect()

def data_unzip():
    # phase 3: unzip data to /cubric/collab/108_QA/SCANNER/raw/XNATID/SESSIONID
    print('\nChecking for XNAT zip files in ' + data_path)
    for ppt, ppt_xn in qcsubj.items():
        #set up temporary dicom directory (scanner level)
        dicom_temp = os.path.join(data_path,ppt,'dicom_temp')
        if not os.path.isdir(dicom_temp):
            os.mkdir(dicom_temp)
            
        # set up temporary nifti directory (scanner level)    
        nifti_temp = os.path.join(data_path, ppt,'nifti')
        if not os.path.isdir(nifti_temp):
            os.mkdir(nifti_temp)        
        # directory of downloaded zip files
        dir_zip = os.path.join(data_path, ppt, 'zip')
        fls = os.listdir(dir_zip)
        for ff in fls:
            if 'XNAT' in ff[0:4] and '.zip' in ff[-4:]:
                exp_id = ff[:-4]
                zip_file = os.path.join(dir_zip,ff)
                # data structure is dir_xnat_sess/dir_scan_sess/scans/
                dir_xnat_sess = os.path.join(dicom_temp,ff[:-4])
                
                # skip if unpacked directory exists
                if os.path.isdir(dir_xnat_sess): 
                    print(dir_xnat_sess + ' exists.  Skipping unzip')
                else:
                    print('Unzipping ' + zip_file + ' to ' + dir_xnat_sess)
                    sb = subprocess.run(['unzip', '-q', '-d', dir_xnat_sess, zip_file])
                
                # get launch dir & output dir for dcm2niix
                dir_scan_sess = os.path.join(dir_xnat_sess,os.listdir(dir_xnat_sess)[0],'scans')
                dir_scan_nifti = os.path.join(nifti_temp, exp_id)
                
                # skip if nifti directory exists
                if os.path.isdir(dir_scan_nifti):
                    print(dir_scan_nifti + ' exists.  Skipping dcm2niix')
                else:   
                    os.mkdir(dir_scan_nifti)
                    print('Running dcm2niix and outputing to  ' + dir_scan_nifti)
                    sb = subprocess.run(['dcm2niix', \
                                     '-f', '%i_%s_%d',\
                                     '-o', dir_scan_nifti,
                                     dir_scan_sess] , \
                                     stdout=subprocess.DEVNULL)
        # tidy up temp directories
        print('Tidying temporary directories' + dicom_temp + ' and ' + dir_zip)
        shutil.rmtree(dicom_temp)
        shutil.rmtree(dir_zip)

#update_download_list()
xnat_download()
data_unzip()