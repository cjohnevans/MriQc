#!/bin/bash
# find summary files from the past day, send to QA channel on MR teams
newfiles=`find /cubric/collab/108_QA -ipath '*glover*fmriqc_latest.txt' -ctime -1 | xargs`
# always send, even if old
newfiles=`find /cubric/collab/108_QA -ipath '*glover*fmriqc_latest.txt'  | xargs`

for ff in $newfiles
do
    name_short=`echo $ff | awk -F ary/ '{print $2}' | awk -F. '{print $1}'`
    name_noext=`echo $ff | awk -F. '{print $1}'`
    name_plot=`echo $name_noext`.png
    scanner=`echo $name_short | awk -F_ '{print $1}'`
    
# to teams
# need to add plots and format title better
#    if [ "$scanner" == "3TW" ]; then  
#        echo "mailx -s $scanner -a $name_plot 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms  < $ff";
#    elif [ "$scanner" == "3TE" ]; then 
#        mailx -s $scanner -a $name_plot 3348574b.cf.onmicrosoft.com@emea.teams.ms  < $ff;
#    elif [ "$scanner" == "7T" ]; then 
#        mailx -s $scanner -a $name_plot d6e057f1.cf.onmicrosoft.com@emea.teams.ms  < $ff;
#    elif [ "$scanner" == "3TM" ]; then 
#        mailx -s $scanner -a $name_plot 76187b6b.cf.onmicrosoft.com@emea.teams.ms  < $ff;
#    fi

# to me
    if [ "$scanner" == "3TW" ]; then  
        mailx -s $scanner -a $name_plot evansj31@cardiff.ac.uk  < $ff;
    elif [ "$scanner" == "3TE" ]; then 
        mailx -s $scanner -a $name_plot evansj31@cardiff.ac.uk  < $ff;
    elif [ "$scanner" == "7T" ]; then 
        mailx -s $scanner -a $name_plot evansj31@cardiff.ac.uk  < $ff;
    elif [ "$scanner" == "3TM" ]; then 
        mailx -s $scanner -a $name_plot evansj31@cardiff.ac.uk  < $ff;
    fi



done

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


