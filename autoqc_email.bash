#!/bin/bash

# FMRIQC
# explicit definition of png files.  use imagemagick's convert to merge the pngs into a single file 
convert -append /cubric/collab/108_QA/QA3TM/fmriqc_glover/summary/3TM_fmriqc_latest.png /cubric/collab/108_QA/QA7T/fmriqc_glover/summary/7T_fmriqc_latest.png /cubric/collab/108_QA/QA3TE/fmriqc_glover/summary/3TE_fmriqc_latest.png  /cubric/collab/108_QA/QA3TW/fmriqc_glover/summary/3TW_fmriqc_latest.png  /cubric/collab/108_QA/fmriqc_latest.png
/home/sapje1/code/python_mrobjects/qc/autoqc_status.py | mailx -r mriqc@cardiff.ac.uk -s "fmriqc (Glover)" -a /cubric/collab/108_QA/fmriqc_latest.png evansj31@cardiff.ac.uk

# QUICKQC
convert -append /cubric/collab/108_QA/QA3TM/quick_SNR_gre3D/summary/3TM_quickqc_latest.png /cubric/collab/108_QA/QA7T/quick_SNR_gre3D/summary/7T_quickqc_latest.png /cubric/collab/108_QA/QA3TE/quick_SNR_gre3D/summary/3TE_quickqc_latest.png  /cubric/collab/108_QA/QA3TW/quick_SNR_gre3D/summary/3TW_quickqc_latest.png  /cubric/collab/108_QA/quickqc_latest.png
/home/sapje1/code/python_mrobjects/qc/autoqc_status.py | mailx -r mriqc@cardiff.ac.uk -s "mriqc (QuickQc)" -a /cubric/collab/108_QA/quickqc_latest.png evansj31@cardiff.ac.uk

spikefiles=`find /cubric/collab/108_QA -ipath '*spike_stats.png' -ctime -1 | xargs`

# loop over spikefiles, and email the new ones
for ss in $spikefiles
do
    name_root=`echo $ss | awk -Fstats '{print $1}'`
    name_stats=$ss
    name_image=$name_root'images.png'
    title=`echo $name_root | awk -F/ '{print $9}'`
    scanner=`echo $title | awk -F_ '{print $1}'`

    if [ "$scanner" == "QA3TW" ]; then  
        echo "spike" | mailx -s $title -a $name_stats -a $name_image 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms;
    elif [ "$scanner" == "QA3TE" ]; then 
        echo "spike" | mailx -s $title -a $name_stats -a $name_image 3348574b.cf.onmicrosoft.com@emea.teams.ms;
    elif [ "$scanner" == "QA7T" ]; then 
        echo "spike" | mailx -s $title -a $name_stats -a $name_image d6e057f1.cf.onmicrosoft.com@emea.teams.ms; 
    elif [ "$scanner" == "QA3TM" ]; then 
        echo "spike" | mailx -s $title -a $name_stats -a $name_image 76187b6b.cf.onmicrosoft.com@emea.teams.ms; 
    fi

done


