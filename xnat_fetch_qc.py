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
sys.path.append('/home/sapje1/code/MriQc')
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

allowed_qc_type = ['fmriqc', 'spikehead', 'spikebody']

def update_downloaded(qc_type):
    """
    update_downloaded()
    ----------------------
    
    Check data downloaded, unzipped and converted to nifti in DATA_PATH/nifti and update 
    file DOWNLOAD_DONE with a set of XNAT experiment IDs.  
    This is to be used by xnat_download() to grab only NEW data from XNAT.
    It should be run before any other calls in this module.
    
    DOWNLOAD_DONE is rewritten on each call (so will update for different data types, e.g.
    when called for spike_head after fmriqc, etc)
    
    Parameters
    ----------
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'
 
    Returns
    -------
    None.

    """
    
    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  update_downloaded() failed")
        return False
    
    # regenerate list of downloaded datasets
    exp_done = []
    
    print('update_downloaded: Checking for downloaded nifti data in ' + data_path)
    for qd in qcsubj.keys():
        qdir = os.path.join(data_path,qd,'nifti',qc_type)
        list1 = os.listdir(qdir)
        exp_done_ppt = []
    
        for l in list1:
            # this checks for the existance of the XNAT_E directory
            if 'XNAT_E' in l and '.zip' not in l:
                # needs to have at least one file
                if len(os.listdir(os.path.join(qdir,l))) > 0:
                    exp_done.append(l)
                    exp_done_ppt.append(l)
        print(qd + ': Found ' + str(len(exp_done_ppt)) + ' downloaded nifti datasets in ' + qdir)
        exp_done_ppt.sort()
        #print(exp_done_ppt)

    #print('Found ' + str(len(exp_done)) + ' total downloaded datasets')
    with open(download_done,'w') as f:
        for line in exp_done:
            f.write(line + '\n')
    print('File ' + download_done + ' updated')

def update_xnat_new(qc_type):
    """
    update_xnat_new()
    ---------------
    
    Check sessions on XNAT matching 108_QA and one of qasubj (QA3TW, ...), and get the xnat experiment
    IDs (XNAT_E...).  Put those that are NOT in DOWNLOAD_DONE into DOWNLOAD_NEW


    Parameters
    ----------
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'
 
    Returns
    -------
    None.

    """
    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  update_xnat_new() failed")
        return False
    if qc_type == 'fmriqc':
        qc_series_name = 'GloverGSQAP'
    if qc_type == 'spikehead':
        qc_series_name = 'EPIspike_head'
    if qc_type == 'spikebody':
        qc_series_name = 'EPIspike_body'
       
 
#   load the list of downloaded and converted sessions
    with open(download_done, 'r') as f:
        tmp = f.readlines()
    exp_downloaded = []
    for line in tmp:
        exp_downloaded.append(line.replace('\n',''))
        
    # if no authentication provided, xnat will look to .netrc file for authentication
    print('\nupdate_xnat_new: Checking for new '+ qc_series_name + ' data on XNAT')
    
    with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
        with open(download_new, 'w') as fnew:
            fnew.write('# 108_QA2023\n') # overwrites old file
        subj = session.projects['108_QA2023'].subjects
        for s in subj:
            for (qc_subj_name, qc_subj_xn) in qcsubj.items():
                if qc_subj_xn in s:
                    subj_id = subj[s].label  # subj_id is QA3TM etc.. 
                    qc_exp = subj[s].experiments
                    # e is all qc experiments for a given subject (=scanner)
                    exp_to_download = []
                    # set up zip path for this scanner (subject)
                    dir_zip = os.path.join(data_path, subj_id, 'zip',qc_type)
                    if not os.path.isdir(dir_zip):
                        os.mkdir(dir_zip)                     
                    
                    for e in qc_exp:
                        # check xnat experiment ids against those previously downloaded
                        if e not in exp_downloaded:
                            if qc_series_name in qc_exp[e].scans:
                                #print('Session ' + e + ' is new.')
                                exp_to_download.append(e)
                            
                    print(subj_id + ': ' + str(len(qc_exp)) + ' sessions on XNAT, ' + \
                          str(len(exp_to_download)) + ' new session(s) are available')
                   
                    with open(download_new, 'a') as fnew:
                        fnew.write('#' + qc_subj_name + '\n')
                        fnew.writelines(line+'\n' for line in exp_to_download)
    session.disconnect()
    print('File ' + download_new + ' updated')


def xnat_download(qc_type, check_all_exams=False):
    """
    xnat_download(check_all=False)
    
    Check through list of XNAT experiments in DOWNLOAD_NEW and download from XNAT.
    Logic of this is updated, following update of XNAT.  Can no longer get all QA
    relevant sequences in a single zip file, need to split up the data in different 
    data downloads.  Download of full exam was failing, but series-by-series download
    works OK.  It does make the directory structure more ungainly, and it makes to 
    download process less efficient, as it needs to be re-run for each data type separately
    
    Parameters
    ----------
    
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'
    
    check_all_exams:  boolean.  
                      If True, check all exams in 108QA project, not just the ones
                      detected as new by update_downloaded().   
    
    """
    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  xnat_download() failed")
        return False  
    if qc_type == 'fmriqc':
        qc_series_name = 'GloverGSQAP'
    if qc_type == 'spikehead':
        qc_series_name = 'EPIspike_head'
    if qc_type == 'spikebody':
        qc_series_name = 'EPIspike_body'
    
    # get list of new experiments to download   
    with open(download_new, 'r') as f:
        tmp = f.readlines()
        exp_new = []
        for line in tmp:
            if 'XNAT_E' in line:
                exp_new.append(line.replace('\n',''))
        print('There are ' + str(len(exp_new)) + ' experiments to check')
                 
    # download new experiments
    with xnat.connect('https://xnat.cubric.cf.ac.uk') as session:
        subj = session.projects['108_QA2023'].subjects
        for s in subj:
            for (qc_subj_name, qc_subj_xn) in qcsubj.items():
                if qc_subj_xn in s:
                    subj_id = subj[s].label  # subj_id is QA3TM etc.. 
                    qc_exp = subj[s].experiments #qc_exp is the qc experiments for this subject
                    # set up zip path for this scanner (subject)
                    dir_zip = os.path.join(data_path, subj_id, 'zip',qc_type)
                    if not os.path.isdir(dir_zip):
                        os.mkdir(dir_zip)                     
                    
                    for ed in qc_exp:   #loop over xnat experiments for this subject
                        if ed in exp_new or check_all_exams == True:  # download if new
                            zip_path = os.path.join(dir_zip, ed + '.zip')
                            print('Checking ' + ed )
                            #print(qc_exp[ed].scans)
                            try:
                                qc_exp[ed].scans[qc_series_name].download(zip_path)
                            except:
                                print('WARNING: No ' + qc_series_name + ' found in ' + ed)
  

def data_unzip(qc_type, unzip=True, remove_invalid_file=False, unzip_all=False):
    """
    data_unzip()
    ------------
    
    Parameters
    ----------
    
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'
   
    unzip          :  Boolean
                      if True performs the extraction, if False, runs a zip file test (unzip -t)
                      without extracting the dicoms.
      
    remove_invalid_file : Boolean
                      if True, remove the .zip file following a failed check or attempted unzip
                      if False, ignore the failed zip file.
       
    unzip_all       : Boolean
       if True, unzip all zip files in the zip directory, irrespective of whether 
       there is already nifti data present (may be useful for exams which have
       multiple data types)
    
    Following on from xnat_download(), check the experiment .zip files in DATA_PATH/SUBJECT_NAME/zip
    then unzip and convert to nifti
    
    Upon conversion to nifti, data are placed in 
    DATA_PATH/SUBJECT_NAME/nifti/QC_TYPE/XNAT_EXPERIMENT_ID
    ( e.g.  /cubric/collab/108_QA/3TE/nifti/QC_TYPE/XNAT_E11630 )
    
    After completion, remove the temporary directories with dicoms (dicom_temp)
    and the zip directory (zip)

    Returns
    -------
    None.

    """
    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  xnat_download() failed")
        return False  
    
    print('\ndata_unzip:  Checking for XNAT zip files in ' + data_path)
    # clean temp directories first
    clean_temp_directories()
    
    for ppt, ppt_xn in qcsubj.items():
        #set up temporary dicom directory (scanner level)
        dir_zip = os.path.join(data_path, ppt, 'zip', qc_type)
        dicom_root = os.path.join(data_path,ppt,'dicom_temp')
        nifti_root = os.path.join(data_path, ppt,'nifti', qc_type)
        if not os.path.isdir(dicom_root):
            os.mkdir(dicom_root)           
        fls = os.listdir(dir_zip)
        for ff in fls:
            if 'XNAT' in ff[0:4] and '.zip' in ff[-4:]:
                exp_id = ff[:-4]
                zip_file = os.path.join(dir_zip,ff)
                # data structure is dicom_exp_dir/dicom_scan_dir/scans/
                dicom_exp_dir = os.path.join(dicom_root,exp_id)
                nifti_exp_dir = os.path.join(nifti_root, exp_id)
                # perform unzipping, but only if unzip_all is True OR
                #  either the dicom or nifti directory exists.
                if unzip_all==True or not os.path.isdir(dicom_exp_dir) or not os.path.isdir(nifti_exp_dir): 
                    if unzip == True:
                        # suppress output - errors are verbose
                        sb = subprocess.run(['unzip', '-q', '-d', dicom_exp_dir, zip_file], \
                                        stdout=subprocess.DEVNULL, \
                                        stderr=subprocess.DEVNULL)
                    else: # check, but don't unzip files
                        sb = subprocess.run(['unzip', '-tq', zip_file], \
                                    stdout=subprocess.DEVNULL, \
                                    stderr=subprocess.DEVNULL)    
                            
                    if sb.returncode == 0:
                        # returns 0 if no errors during unzip 
                        print(ff, 'OK        - Unzipping ' + zip_file + ' succeeded')
                    else:
                        if remove_invalid_file:
                            os.remove(zip_file)
                            print(ff, '!!! ERROR - Unzipping ' + zip_file + ' failed. File removed.')
                        else:
                            print(ff, '!!! ERROR - Unzipping ' + zip_file + ' failed. File not removed')
                else:
                    print(ff, 'Skipping  - previous unzip or nifti exists in ', ppt)

def empty_nifti_dir(nifti_dir, remove=False):
    """
    check for null directories in nifti_dir
    return a list of empty directories - XNAT download failures
    """
    dirs = os.listdir(nifti_dir)
    
    empty_dirs = []
    
    for dd in dirs:
        ff = os.listdir(os.path.join(nifti_dir,dd))
        if len(ff) == 0:
            empty_dirs.append(dd)
            if remove == True:
                os.rmdir(os.path.join(nifti_dir,dd))
    print(nifti_dir + ' has ' + str(len(dirs)) + ' directories, and ' + str(len(empty_dirs)) + ' empty directories')

    if len(empty_dirs) > 0 and remove == True:
        print("Empty directories removed")

def clean_temp_directories():
    """
    clean_temp_directories()
    
    Clean up temporary directories to try to prevent reanalysis of
    data when switching between branches

    Returns
    -------
    None.

    """
    for qd in qcsubj.keys():
        d=os.path.join(data_path,qd,'dicom_temp')
        print('removing ' + d)
        subprocess.run(['rm', '-r', d])
        # recreate now, in case not running analysis in order 
        subprocess.run(['mkdir', d])

def nifti_convert(qc_type):
    '''
    Convert files in dicom_temp dir to nifti.  
    Check for existing nii files in relevant nifti directory (from previous conversion),
    and skip if nii file already exists for this exam.
    
    
    Keep dicom_temp clean - delete after an unpacking attempt.
    
    Parameters
    ----------
    
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'

    Returns
    -------
    None.

    '''

    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  xnat_download() failed")
        return False  
     
    print('\nnifti_convert: Checking for dicom data in ' + data_path)
    
    for ppt, ppt_xn in qcsubj.items():
        #set up temporary dicom directory (scanner level)
        dicom_root = os.path.join(data_path, ppt,'dicom_temp')
        nifti_root = os.path.join(data_path, ppt,'nifti', qc_type)
        empty_nifti_dir(nifti_root, remove=True) # clean the nifti dir before starting

        fls = os.listdir(dicom_root)
        for ff in fls:
            if 'XNAT' in ff[0:4]:
                exp_id = ff
                dicom_exp_dir = os.path.join(dicom_root,exp_id)
                dicom_scan_dir = os.path.join(dicom_exp_dir,os.listdir(dicom_exp_dir)[0],'scans')
                nifti_exp_dir = os.path.join(nifti_root, exp_id)

                # check if niftis exist - sufficient as we've cleaned up any empty dirs.
                if os.path.isdir(nifti_exp_dir): 
                    print(ff, 'Skipping  - previous unzip or nifti exists in ', ppt)
                else:
                    os.mkdir(nifti_exp_dir)
                    #print(ff, 'Running dcm2niix and outputing to  ' + nifti_exp_dir)
                    sb = subprocess.run(['dcm2niix', \
                            '-f', '%i_%s_%d',\
                            '-o', nifti_exp_dir,
                            dicom_scan_dir] , \
                            stdout=subprocess.DEVNULL)
                    if sb.returncode == 0:
                        # returns 0 if no errors during unzip 
                        print(ff, 'OK        - Convert from ' + dicom_exp_dir + ' succeeded')                        
                    else:
                        print(ff, '!!! ERROR - Convert ' + dicom_exp_dir + ' failed.  Error=' + sb.returncode)
        # clear d        
             
def proc_qc(qc_type, analyse_all=False):
    """
    proc_fmriqc(analyse_all=False):

    Principle:
       get the nifti names in names in scanner/nifti directories
       filter the fmri runs
       get the names of the report directories in the scanner/summary directories
       only analyse the niftis where there is no existing report        

    Parameters
    ----------
    qc_type        :  str
                      Must be one of allowed_qc_type (defined above), e.g. 'fmriqc', 'spikehead'
    
    analyse_all    :  BOOL, optional
                      Reanalyse all data. If True then reanalyse all data, irrespective
                      of whether output already exists.

    Returns
    -------
    None.

    """
    
    if qc_type not in allowed_qc_type:
        print(qc_type + " is not an allowed qc type.  xnat_download() failed")
        return False  
    if qc_type == 'fmriqc':
        qc_series_name = 'GloverGSQAP'
    if qc_type == 'spikehead':
        qc_series_name = 'EPIspike_head'
    if qc_type == 'spikebody':
        qc_series_name = 'EPIspike_body'
    
    
    # set up paths
    data_path = '/cubric/collab/108_QA'
    scanners = ['QA7T', 'QA3TW', 'QA3TE', 'QA3TM']
    nifti_path = []

    for s in scanners:
        print('Checking for new ' + qc_type + ' data in '+ s)
        nifti_path = os.path.join(data_path, s, 'nifti', qc_type)
        exam_list = os.listdir(nifti_path)
        for e in exam_list:
            # find fmriqc data types
            series_list = os.listdir(os.path.join(nifti_path, e))
            for se in series_list:
                se_no_ext = se.split('.')[0] # filename, no extension
                scan_date = se_no_ext[0:17]
                if qc_series_name in se and '.nii' in se:
                    rep_path = os.path.join(data_path, s, qc_type, 'proc',se_no_ext)
                    try:
                        os.listdir(rep_path)
                    except:
                        print(se_no_ext + ' is new... Analysing...')
                        if qc_type == 'fmriqc':
                            mriqc.FmriQc(os.path.join(nifti_path, e, se) \
                                  , report_path = rep_path)
                        if qc_type == 'spikehead' or qc_type == 'spikebody':
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
