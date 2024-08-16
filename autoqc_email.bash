#!/bin/bash

# ------------- fmriqc Summary ---------------------
# use imagemagick's convert to merge the pngs into a single file 
# remove 3TW from the fmriqc summary
convert -append /cubric/collab/108_QA/QA7T/fmriqc/summary/7T_fmriqc_latest.png \
  /cubric/collab/108_QA/QA3TM/fmriqc/summary/3TM_fmriqc_latest.png \
  /cubric/collab/108_QA/QA3TE/fmriqc/summary/3TE_fmriqc_latest.png  \
  /cubric/collab/108_QA/fmriqc_latest.png

/home/sapje1/code/MriQc/autoqc_status.py | \
  mailx -r mriqc@cardiff.ac.uk -s "fmriqc Summary" -a /cubric/collab/108_QA/fmriqc_latest.png \
  5ed88dc2.cf.onmicrosoft.com@emea.teams.ms


# ------------- fmriqc recent ---------------------
recent=`find /cubric/collab/108_QA -ipath '*fmriqc*proc*GloverGSQAP' -ctime -1`
for ff in $recent
do
    # save as jpeg to prevent replicating images in multiple runs.
    combine=`echo $ff | awk -F/ '{print $8}' | awk -F- '{print $1"-"$2"_fmriqc.jpg"}'`
    outfile=${ff}'/'${combine}
    convert -append $recent/*.png $outfile
#    convert -append $recent/SFNR.png $recent/pixel_histogram.png \
#      $recent/slice_time.png $recent/drift_correct.png $outfile
# for some reason, this one doesn't like the reply-to field??
    title="fmriqc Recent"$combine
    echo "Recent fmriqc report from "$ff | \
     mailx -s "fmriqc Recent" -a $outfile 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms
#    echo "Recent fmriqc report" | \
#     mailx -r mriqc@cardiff.ac.uk -s "fmriqc Recent" -a $outfile 5ed88dc2.cf.onmicrosoft.com@emea.teams.ms

done

# ------------- spike recent ---------------------
#spikefiles=`find /cubric/collab/108_QA -ipath '*spike_stats.png' -ctime -1 | xargs`
spikefiles=`find /cubric/collab/108_QA -ipath '*spike_stats.png' -mmin -240 | xargs`

# loop over spikefiles, and email the new ones
for ss in $spikefiles
do
    name_root=`echo $ss | awk -Fstats '{print $1}'`
    name_stats=$ss
    name_image=$name_root'images.png'
    title=`echo $name_root | awk -F/ '{print $9}'`
    echo "spike" | mailx -s $title -a $name_stats -a $name_image \
     5ed88dc2.cf.onmicrosoft.com@emea.teams.ms;
done


