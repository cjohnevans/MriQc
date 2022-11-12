#!/bin/bash

qcdir='/home/sapje1/data_sapje1/QC/RoutineQC/'
transferdir="${qcdir}summary/fortransfer"


for txtfile in `ls ${transferdir}/QA*.txt`
  do basefname=`echo ${txtfile} | awk -Ftransfer/ '{print $2}' | awk -F. '{print $1}'`
  tail -n 1 $txtfile | awk -F, '{print "Scanner\t\t"$2 "\nDate\t\t" $3 "\nTime\t\t" $4 "\nSNR\t\t" $5 "\nSNRvol\t\t" $6 "\nSFNRroi\t\t" $7 "\nSFNRroi_nodrift\t" $8 "\nLinearDrift\t" $11 "\nSliceMaxSigChg\t" $12 "\nDriftArtefact\t" $13}' > "$transferdir/${basefname}.dat"
done 

# Send single session QA report to Quality Control - Scan Reports channel
# Quality Control - Scan Reports - MR Team <d6e057f1.cf.onmicrosoft.com@emea.teams.ms>
for datfile in `ls ${transferdir}/QA*.dat`
  do basefname=`echo ${datfile} | awk -Ftransfer/ '{print $2}' | awk -F. '{print $1}'`
  pdfname="${transferdir}/${basefname}.pdf"
  subj=`echo $basefname | awk -F- '{print $1 "___" $3 "___" $4}'`
  cat ${datfile} | mailx -s $subj -a $pdfname d6e057f1.cf.onmicrosoft.com@emea.teams.ms
done

# Send Shewhart reports to Quality Control channel
# Quality Control - MR Team <5ed88dc2.cf.onmicrosoft.com@emea.teams.ms>
echo "3T West" | mailx -s "3T West Shewhart QC" -a ${transferdir}/Shewhart3TW.pdf 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
echo "3T East" | mailx -s "3T East Shewhart QC" -a ${transferdir}/Shewhart3TE.pdf 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
echo "Connectom" | mailx -s "Connectom Shewhart QC " -a ${transferdir}/ShewhartConnectom.pdf 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
echo "7T" | mailx -s "7T Shewhart QC" -a ${transferdir}/Shewhart7T.pdf 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
