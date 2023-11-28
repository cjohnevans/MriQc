#!/bin/bash

newfiles=`find /cubric/collab/108_QA -ipath '*glover*fmriqc_latest.txt' -ctime -1 | xargs`
echo $newfiles

for ff in $newfiles
do
    echo $ff;
# need to add plots and format tiltle better
    mailx -s $ff c.john.evans@gmail.com < $ff;
done

