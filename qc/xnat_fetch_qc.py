#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xnat_fetch_qc.py

Collect qc data from xnat project 108 using xnatpy
Created on Wed Aug 16 14:37:20 2023

Set of functions to handle get data from XNAT and launch processing, with approiate directory structure


@author: sapje1
"""

import xnat
import subprocess
import os, shutil, sys
sys.path.append('/home/sapje1/code/python_mrobjects/qc')
import mriqc

data_path = '/cubric/collab/108_QA'
# download_done is a list of xnat experiments which will be IGNORED by xnat_update_new()
download_done = os.path.join(data_path, 'xnat_download_done.txt')
# download_new is the file created by xnat_update_new which will be downloaded by xnat_download()
download_new =  os.path.join(data_path, 'xnat_new.txt')

qcsubj = { 'QA7T' : 'XNAT_S06014',
           'QA3TM' : 'XNAT_S06051',
           'QA3TE' : 'XNAT_S06097',
           'QA3TW' : 'XNAT_S06091'
        }

def update_downloaded():
    """
    update_downloaded()
    ----------------------
    
    Check data downloaded, unzipped and converted to nifti in DATA_PATH/nifti  and update 
    file DOWNLOAD_DONE with a set
    of XNAT experiment IDs.  This is to be used by xnat_download() to grab only NEW data from XNAT.
    It should be run before any other calls in this module.

    Returns
    -------
    None.

    """
    # regenerate list of downloaded datasets
    exp_done = []
    
    print('Checking for downloaded data in ' + data_path)
    for qd in qcsubj.keys():
        qdir = os.path.join(data_path,qd,'nifti')
        list1 = os.listdir(qdir)
        exp_done_ppt = []
    
        for l in list1:
            # this checks for the existance of the XNAT_E directory
            if 'XNAT_E' in l and '.zip' not in l:
                # needs to have at least one file
                if len(os.listdir(os.path.join(qdir,l))) > 0:
                    exp_done.append(l)
                    exp_done_ppt.append(l)
        print('Found ' + str(len(exp_done_ppt)) + ' downloaded datasets for ' + qd)
        exp_done_ppt.sort()
        print(exp_done_ppt)

    print('Found ' + str(len(exp_done)) + ' total downloaded datasets')
    with open(download_done,'w') as f:
        for line in exp_done:
            f.write(line + '\n')

def update_xnat_new(update_local=True):
    """
    update_xnat_new()
    ---------------
    
    Check sessions on XNAT matching 108_QA and one of qasubj (QA3TW, ...), and get the xnat experiment
    IDs (XNAT_E...).  Put those that are NOT in DOWNLOAD_DONE into DOWNLOAD_NEW


    Parameters
    ----------
    None

    Returns
    -------
    None.

    """

#   update the list of downloaded and converted sessions
    if update_local:
        update_downloaded()

    with open(download_done, 'r') as f:
        tmp = f.readlines()
    exp_downloaded = []
    for line in tmp:
        exp_downloaded.append(line.replace('\n',''))
        
    # if no authentication provided, xnat will look to .netrc file for authentication
    print('\nChecking for new data on XNAT')
    
    with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
        with open(download_new, 'w') as fnew:
            fnew.write('# 108_QA2023\n') # overwrites old file
        subj = session.projects['108_QA2023'].subjects
        for s in subj:
            for (qc_subj_name, qc_subj_xn) in qcsubj.items():
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
                            
                    print(str(len(exp_to_download)) + ' new session(s) are available\n')
                    with open(download_new, 'a') as fnew:
                        fnew.write('#' + qc_subj_name + '\n')
                        fnew.writelines(line+'\n' for line in exp_to_download)
    session.disconnect()

def xnat_download():
    """
    xnat_download()
    Check through list of XNAT experiments in DOWNLOAD_NEW and download from XNAT. 
    """

    # get list of new experiments to download   
    with open(download_new, 'r') as f:
        tmp = f.readlines()
        exp_new = []
        for line in tmp:
            if 'XNAT_E' in line:
                exp_new.append(line.replace('\n',''))

        print('There are ' + str(len(exp_new)) + ' experiments to download')

                 
    # download new experiments
    with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
        subj = session.projects['108_QA2023'].subjects
        for s in subj:
            for (qc_subj_name, qc_subj_xn) in qcsubj.items():
                if qc_subj_xn in s:
                    subj_id = subj[s].label  # subj_id is QA3TM etc.. 
                    qc_exp = subj[s].experiments #qc_exp is the qc experiments for this subject

                    # set up zip path for this scanner (subject)
                    dir_zip = os.path.join(data_path, subj_id, 'zip')
                    if not os.path.isdir(dir_zip):
                        os.mkdir(dir_zip)                     
                    
                    for ed in qc_exp:   #loop over xnat experiments for this subject
                        if ed in exp_new:  # download if new
                            zip_path = os.path.join(dir_zip, ed + '.zip')
                            print('Downloading ' + ed + ' to ' + zip_path)
                            try:
                                qc_exp[ed].download(zip_path)
                            except:
                                print('!!!WARNING: Download of ' + ed + ' failed')
 
    
def empty_nifti_dir(nifti_dir):
    """
    check for null directories in nifti_dir
    return a list of empty directories - XNAT download failures
    """
    dirs = os.listdir(nifti_dir)
    print(nifti_dir)
    print(dirs)
    
    empty_dirs = []
    
    for dd in dirs:
        ff = os.listdir(os.path.join(nifti_dir,dd))
        if len(ff) == 0:
            empty_dirs.append(dd)
            os.rmdir(os.path.join(nifti_dir,dd))
    print(empty_dirs)
    
    

def data_unzip():
    """
    data_unzip()
    ------------
    
    Following on from xnat_download(), check the experiment .zip files in DATA_PATH/SUBJECT_NAME/zip
    then unzip and convert to nifti
    
    Upon conversion to nifti, data are placed in 
    DATA_PATH/SUBJECT_NAME/nifti/XNAT_EXPERIMENT_ID
    ( e.g.  /cubric/collab/108_QA/3TE/nifti/XNAT_E11630 )
    
    After completion, remove the temporary directories with dicoms (dicom_temp)
    and the zip directory (zip)

    Returns
    -------
    None.

    """
    print('\nChecking for XNAT zip files in ' + data_path)
    for ppt, ppt_xn in qcsubj.items():
        #set up temporary dicom directory (scanner level)
        dicom_root = os.path.join(data_path,ppt,'dicom_temp')
        if not os.path.isdir(dicom_root):
            os.mkdir(dicom_root)
            
        # set up nifti directory (scanner level)    
        nifti_root = os.path.join(data_path, ppt,'nifti')
        if not os.path.isdir(nifti_root):
            os.mkdir(nifti_root)        
        # directory of downloaded zip files
        dir_zip = os.path.join(data_path, ppt, 'zip')
        fls = os.listdir(dir_zip)
        for ff in fls:
            if 'XNAT' in ff[0:4] and '.zip' in ff[-4:]:
                exp_id = ff[:-4]
                zip_file = os.path.join(dir_zip,ff)
                # data structure is dicom_exp_dir/dir_scan_sess/scans/
                dicom_exp_dir = os.path.join(dicom_root,exp_id)
                nifti_exp_dir = os.path.join(nifti_root, exp_id)
                
                # skip if unpacked directory exists
                if os.path.isdir(dicom_exp_dir) or os.path.isdir(nifti_exp_dir): 
                    print(ff, 'Skipping  - previous unzip or nifti exists in ', ppt)
                else:
                    try:
                        # suppress output - errors are verbose
                        sb = subprocess.run(['unzip', '-q', '-d', dicom_exp_dir, zip_file], \
                                            stdout=subprocess.DEVNULL, \
                                            stderr=subprocess.DEVNULL)
                        os.listdir(dicom_exp_dir)
                    except:
                        print(ff, '!!! ERROR - Unzipping ' + zip_file + ' failed')


def nifti_convert():
    '''
    convert files in dicom_temp dir to nifti.  
    work in progress

    Returns
    -------
    None.

    '''
    
    for ppt, ppt_xn in qcsubj.items():
        #set up temporary dicom directory (scanner level)
        dicom_root = os.path.join(data_path, ppt,'dicom_temp')
        nifti_root = os.path.join(data_path, ppt,'nifti')
        fls = os.listdir(dicom_root)

        # need equivalent of this..
        for ff in fls:
            if 'XNAT' in ff[0:4]:

                exp_id = ff

                dicom_exp_dir = os.path.join(dicom_root,exp_id)
                dicom_scan_dir = os.path.join(dicom_exp_dir,os.listdir(dicom_exp_dir[0],'scans'))
                nifti_exp_dir = os.path.join(nifti_root, exp_id)

                # check if niftis exist
                if os.path.isdir(nifti_exp_dir): 
                    print(ff, 'Skipping  - previous unzip or nifti exists in ', ppt)
                else:
                    os.mkdir(nifti_exp_dir)
                print('Running dcm2niix and outputing to  ' + nifti_exp_dir)
                sb = subprocess.run(['dcm2niix', \
                            '-f', '%i_%s_%d',\
                            '-o', nifti_exp_dir,
                            dir_scan_sess] , \
                            stdout=subprocess.DEVNULL)

        
        
def proc_qc(analyse_all=False):
    """
    proc_fmriqc(analyse_all=False):

    Principle:
       get the nifti names in names in scanner/nifti directories
       filter the fmri runs
       get the names of the report directories in the scanner/summary directories
       only analyse the niftis where there is no existing report        


    Parameters
    ----------
    analyse_all : BOOL, optional
        Reanalyse all data. The default is False.  True will reanalyse all data, irrespective
        of whether output already exists.

    Returns
    -------
    None.

    """
    
    # set up paths
    data_path = '/cubric/collab/108_QA'
    scanners = ['QA7T', 'QA3TW', 'QA3TE', 'QA3TM']
    nifti_path = []
    report_path = []

    for s in scanners:
        print('Checking for new nifti files in '+ s)
        nifti_path = os.path.join(data_path, s, 'nifti')
#        report_path = os.path.join(data_path, s, 'fmriqc_glover/proc')
        exam_list = os.listdir(os.path.join(data_path, s, 'nifti'))
        for e in exam_list:
            # find fmriqc data types
            series_list = os.listdir(os.path.join(nifti_path, e))
            for se in series_list:
                se_no_ext = se.split('.')[0] # filename, no extension
                scan_date = se_no_ext[0:8]
                if 'GloverGSQAP.nii' in se:
                    rep_path = os.path.join(data_path, s, 'fmriqc_glover/proc',se_no_ext)
                    try:
                        os.listdir(rep_path)
                    except:
                        print(se_no_ext + ' is new... Analysing...')
                        fmri_qc = mriqc.FmriQc(os.path.join(nifti_path, e, se) \
                                  , report_path = rep_path)
                if 'Warmingup.nii' in se or 'WarmingUp.nii' in se:
                    rep_path = os.path.join(data_path, s, 'fmriqc_warmup/proc', se_no_ext)
                    try:
                        os.listdir(rep_path)
                    except:
                        print(se_no_ext + ' is new... Analysing...')
                        fmri_qc = mriqc.FmriQc(os.path.join(nifti_path, e, se) \
                                  , report_path = rep_path)
                if 'QC_MB_GRE_EPI_FA15_Tx220V.nii' in se:
                    rep_path = os.path.join(data_path, s, 'fmriqc_MB/proc', se_no_ext)
                    try:
                        os.listdir(rep_path)
                    except:
                        print(se_no_ext + ' is new... Analysing...')
                        fmri_qc = mriqc.FmriQc(os.path.join(nifti_path, e, se) \
                                  , report_path = rep_path)
                if 'spike' in se:
                    if '.nii' in se:
                        rep_path = os.path.join(data_path, s, 'spike/proc',se_no_ext)
                        try:
                            os.listdir(rep_path)
                        except:
                            print(se_no_ext + ' is new... Analysing...')
                            plttitle= s + '_' + scan_date
                            mriqc.SpikeQc(os.path.join(nifti_path,e,se) \
                                    , report_path=rep_path).spike_check(plttitle)
    
def check_qc():
    '''
    check the status of the QC directory
    '''
    scanners = ['QA7T', 'QA3TW', 'QA3TE', 'QA3TM']
    for s in scanners:
        print('Summary of ' + s)
        nifti_sess = os.listdir(os.path.join(data_path, s, 'nifti'))
        print('Number of sessions with niftis:'  + str(len(nifti_sess)))


def fetch_qc(download=True, unzip=True, proc_fmri=True):
    update_xnat_new()
    if download:
        xnat_download()
    if unzip:
        data_unzip()
    if proc_fmri:
        proc_qc()
