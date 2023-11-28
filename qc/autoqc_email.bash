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

