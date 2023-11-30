#!/bin/bash
# find summary files from the past day, send to QA channel on MR teams
newfiles=`find /cubric/collab/108_QA -ipath '*glover*fmriqc_latest.txt' -ctime -1 | xargs`

for ff in $newfiles
do
    name_short=`echo $ff | awk -F ary/ '{print $2}' | awk -F. '{print $1}'`
    name_noext=`echo $ff | awk -F. '{print $1}'`
    name_plot=`echo $name_noext`.png
# need to add plots and format tiltle better
#    mailx -s $name_short -a $name_plot c.john.evans@gmail.com < $ff;
    mailx -s $name_short -a $name_plot d6e057f1.cf.onmicrosoft.com@emea.teams.ms  < $ff;
done

spikefiles=`find /cubric/collab/108_QA -ipath '*spike_stats.png' -ctime -1 | xargs`

# loop over spikefiles, and email the new ones
for ss in $spikefiles
do
    name_root=`echo $ss | awk -Fstats '{print $1}'`
    name_stats=$ss
    name_image=$name_root'images.png'
    title=`echo $name_root | awk -F/ '{print $9}'`
    echo "spike" | mailx -s $title -a $name_stats -a $name_image d6e057f1.cf.onmicrosoft.com@emea.teams.ms;
done


