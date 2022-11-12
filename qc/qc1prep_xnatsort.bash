#!/bin/bash
# qc0prep_xnatsort
#  1) Check contents of $XNATZIPDIR directory for XNAT downloaded images (*.zip)
#  2) Unzip all .zip files in here
#  3) Convert to nifti (dcm2niix)
#  4) move GloverEPI, spike, phase stability and ghosting niftis into separate dirs, leave 
#     others in 'uncategorised' 

QCDIR="/home/sapje1/data_sapje1/QC/RoutineQC"
XNATZIPDIR=$QCDIR"/xnatzip_input"
XNATUNZIPDIR=$QCDIR"/proc_temp/xnat_unzip"  # unzipped xnat data
PROCDIR=$QCDIR"/proc_temp/nii_proc"         # data to be used
UNCATDIR=$QCDIR"/proc_temp/nii_ignore"   # uncategorised - data to be ignored

cd $XNATZIPDIR

# check for unzipped data from XNAT first
for dir in `ls -d */`; do 
  echo "Moving " $dir
  mv $dir $XNATUNZIPDIR
done

for file in `ls *.zip`; do
  scandir=`echo $file | sed 's/.zip//'`
  echo "Unzipping " $file
  unzip -q $file -d $XNATUNZIPDIR/$scandir
done


/cubric/software/bin/dcm2niix -z n -v n -f %n-%d-%i -o $UNCATDIR $XNATUNZIPDIR
for niigzfilepath in $UNCATDIR/*Warming* ; do mv $niigzfilepath $PROCDIR; done

