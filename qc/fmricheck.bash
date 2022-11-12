#!/bin/bash
# input : fmricheck INPUT_NII_GZ
#         need the input dataset to be named in first argument (nii.gz)
#         don't include the .nii.gz in INPUT_NII_GZ
# output: output_CoV.nii.gz (map of coefficient of variation)
#         output_SFNR.nii.gz (map of SFNR)
#         output_mean_CoV.txt (file with mean CoV for non-zero voxels)
#         output_mean_SFNR.txt (file with mean SFNR for non-zero voxels)

mcflirt -in $1.nii.gz -out  $1_mc -stats
bet $1_mc.nii.gz $1_mask -m
fslmaths $1_mc.nii.gz -mul $1_mask_mask $1_mc_brain
fslmaths $1_mc_brain.nii.gz -Tstd $1_mc_brain_Tstd
fslmaths $1_mc_brain.nii.gz -Tmean $1_mc_brain_Tmean
fslmaths $1_mc_brain_Tstd -div $1_mc_brain_Tmean -mul 100 $1_output_CoV
fslmaths $1_mc_brain_Tmean -div $1_mc_brain_Tstd.nii.gz $1_output_SFNR
fslstats $1_output_CoV.nii.gz -M > $1_output_mean_CoV.txt
fslstats $1_output_SFNR.nii.gz -M > $1_output_mean_SFNR.txt

