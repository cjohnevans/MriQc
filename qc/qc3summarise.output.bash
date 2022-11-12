#!/bin/bash
# plan
#  1) generate summary text files for four scanners
#  ?2) run an automated matlab script to do the autoreporting?


SUMMARYDIR='/home/sapje1/data_sapje1/QC/RoutineQC/summary'
PROCOUTDIR='/home/sapje1/data_sapje1/QC/RoutineQC/proc_output'

#get the header
headertxt=`cat ${PROCOUTDIR}/*.txt | head -n 1`

# generate 3TE summary file
echo $headertxt > $SUMMARYDIR/3TEQC_epi.txt
for epifile in $PROCOUTDIR/QA3TE*.txt;
  do tail -n 1 $epifile >> $SUMMARYDIR/3TEQC_epi.txt
done

# generate 3TW summary file
echo $headertxt > $SUMMARYDIR/3TWQC_epi.txt
for epifile in $PROCOUTDIR/QA3TW*.txt;
  do tail -n 1 $epifile >> $SUMMARYDIR/3TWQC_epi.txt
done

# generate Connectom summary file
echo $headertxt > $SUMMARYDIR/ConnectomQC_epi.txt
for epifile in $PROCOUTDIR/QA3TM*.txt;
  do tail -n 1 $epifile >> $SUMMARYDIR/ConnectomQC_epi.txt
done

# generate 7T summary file
echo $headertxt > $SUMMARYDIR/7TQC_epi.txt
for epifile in $PROCOUTDIR/QA7T*.txt;
  do tail -n 1 $epifile >> $SUMMARYDIR/7TQC_epi.txt
done

