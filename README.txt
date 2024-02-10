CJE 22/4/2020
Quick Start:
1) place xnat data (or links to zip in archive) in either 
    a) the xnat zipped data in xnatzip_input dir, or
    b) the unzipped xnat data in proc_temp/xnat_unzip
2) run routineqc.bash, this will...
      run qc1prep_xnatsort.bash
      run qc2analyse.epi.bash
      run qc3summarise.output.bash
      run qc3summarise.cleanup.bash
3) from matlab, run qc3_summarise_plotQC.m 

Data Structure:
proc_nii: temp directory with niftis to be processed.
summary: output summary data to be stored (PDFs, summary text files)
xnatzip_input : temp directory where the .tar.gz files to be processed should be placed.
zip_archive: storage of XNAT .tar.gz files

/cubric/data/sapje1/QA/RoutineQC
xnatzip_input: temp directory where the .tar.gz files to be processed should be placed.
proc_output: output summary data from individual datasets (PDFs, summary text files)
summary: summary data across multiple datasets

/cubric/scratch/sapje1/RoutineQC
proc_temp: temp directory used during processing
xnatzip_archive: archive of QC zip data from xnat.


The overall pipeline for RoutineQC is as follows;

----- QC PHASE 0: ACQUIRE/TRANSFER -----
- Data acquired and transferred to XNAT.

- QC data downloaded from XNAT to a 'watched' directory 
  XNATZIPDIR=$QCDIR"/xnat_zip"
This directory is defined in PrepXnatQC.bash
Currently this needs to be a manual process, possibly run by the operators 
after acquiring the QC data, but in future may be something which could be 
automated by XNAT.

----- QC PHASE 1 (qc1) : PREP -----
Scripts;
qc1prep_xnatsort.bash

- script checks XNATZIPDIR directory for any data, and if 
present, unzip and convert to nifti.  Data are required to be in XNAT format,
as downloaded by 'Manage Files'.
Currently this needs to be run manually, but could set up a cron job to 
achieve this.

----- QC PHASE 2 (qc2) : ANALYSE -----
Scripts;
qc2analyse.epi.bash
qc2analyse.spike.bash

- qc2analyse.epi.bash: 
from Matlab run epiQA_Siemens.  This generates a PDF output.  Results in 
matlab are a struct with the format

res = 

                    SNR: 579.1656
                SNR_vol: 1.0993e+03
               SFNR_roi: 1.3237e+03
       SFNR_roi_nodrift: 1.3297e+03
       SFNR_vol_nodrift: 1.2613e+03
    semSFNR_vol_nodrift: 2.1860
            LinearDrift: -0.0221
      SliceMaxSigChange: 0.1193
          DriftArtefact: 0.6483

Need to store these results in a results file for each output.

-- QC PHASE 3 (qc3) : SUMMARISE --
Scripts to generate long term plots of changes over time.  Not yet complete
qc3summarise.output.bash:  look through all epiqc summaru .txt files and compile these into a single output
file per scanner in the summary directory.  This will be used by the qc3_summarise_plotQC.m matlab script.

qc3summarise.cleanup.bash: tidy up directories, ready for next run.

qc3_summarise_plotQC.m: Run the Shewhart analysis on the summary data produced by qc3summarise.output.bash, and send the resulting plots to the MS Teams channel MR - Quality Control.


