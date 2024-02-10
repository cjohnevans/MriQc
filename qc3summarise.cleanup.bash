#!/bin/bash
#qc3summarise.cleanup

QCDIR="/home/sapje1/data_sapje1/QC/RoutineQC"
XNATZIPDIR=$QCDIR"/xnatzip_input"
XNATUNZIPDIR=$QCDIR"/proc_temp/xnat_unzip"  # unzipped xnat data
PROCDIR=$QCDIR"/proc_temp/nii_proc"         # data to be used
UNCATDIR=$QCDIR"/proc_temp/nii_ignore"   # uncategorised - data to be ignored

rm $UNCATDIR/* 2> /dev/null
rm $PROCDIR/* 2> /dev/null
rm $XNATZIPDIR/* 2> /dev/null
rm -r $XNATUNZIPDIR/* 2> /dev/null
rm -r $QCDIR/summary/fortransfer/*
