#!/bin/bash
datdir="/home/sapje1/data_sapje1/2022/fmri_param_check/nifti/"
echo $datdir
rm $datdir/tstd.txt

for file in $datdir/*a???.nii.gz;
  do filenoext=`echo $file |awk -F.nii '{print $1}'`
  echo ${filenoext};
  echo bet ${filenoext} ${filenoext}_bet -m;
  bet ${filenoext} ${filenoext}_bet -m;
  mcflirt -in ${filenoext}_bet
  echo fslmaths ${filenoext}_bet_mcf -Tmean ${filenoext}_tmean;  
  fslmaths ${filenoext}_bet_mcf -Tmean ${filenoext}_tmean;  
  echo fslmaths ${filenoext}_bet_mcf -Tstd ${filenoext}_tstd;
  fslmaths ${filenoext}_bet_mcf -Tstd ${filenoext}_tstd;
  echo fslmaths ${filenoext}_tmean -div ${filenoext}_tstd tmp;
  fslmaths ${filenoext}_tmean -div ${filenoext}_tstd tmp;
  echo fslmaths tmp -mul ${filenoext}_bet_mask ${filenoext}_tsnr;
  fslmaths tmp -mul ${filenoext}_bet_mask ${filenoext}_tsnr;
  fslstats ${filenoext}_tsnr -M
  rm tmp.nii.gz 
done

for file in $datdir/*a???_tsnr.nii.gz;
  do echo ${file} >> ${datdir}/tsnr.txt;  
  fslstats ${file} -M >> ${datdir}/tsnr.txt;
done
