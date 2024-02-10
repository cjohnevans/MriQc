#!/bin/bash 
# wrapper script for QC

QCDIR="/home/sapje1/data_sapje1/QC/RoutineQC"
XNATZIPDIR=$QCDIR"/xnatzip_input"
routineqcpath=/home/sapje1/code/QA/RoutineQC

echo "------------------------------------------------------------------------"
echo "  routineqc.bash (no arguments)"
echo "      check xnatzip_input dir, process any input data  then run summary"
echo "      scripts.  If no zip files present, summary scripts are still run"
echo
echo "  routineqc.bash cron"
echo "      for crontab.  If no input data, do nothing."
echo "------------------------------------------------------------------------"

newdata=false
forcesummary=true
if [ "$1" == "cron" ]; then
   forcesummary=false
fi

if [ "$(ls -A ${XNATZIPDIR})" ]; then 
   newdata=true
   ${routineqcpath}/qc1prep_xnatsort.bash
   ${routineqcpath}/qc2analyse_runepiQA.bash
fi

echo "newdata "$newdata
echo "forcesummary "$forcesummary

if [ "$newdata" == true ] || [ "$forcesummary" == true ] ; then
   ${routineqcpath}/qc3summarise.output.bash
   ${routineqcpath}/qc3summarise.plotQC.bash
   ${routineqcpath}/qc3summarise.email.bash
   ${routineqcpath}/qc3summarise.cleanup.bash
fi

