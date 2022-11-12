#!/bin/bash
# input : fmricheck_phantom INPUT_NII_GZ
#         need the input dataset to be named in first argument (nii.gz)
#         don't include the .nii.gz in INPUT_NII_GZ
# output: output_CoV.nii.gz (map of coefficient of variation)
#         output_SFNR.nii.gz (map of SFNR)
#         output_mean_CoV.txt (file with mean CoV for non-zero voxels)
#         output_mean_SFNR.txt (file with mean SFNR for non-zero voxels)

fslmaths $1.nii.gz -Tstd $1_Tstd
fslmaths $1.nii.gz -Tmean $1_Tmean
fslmaths $1_Tstd -div $1_Tmean -mul 100 $1_output_CoV
fslmaths $1_Tmean -div $1_Tstd.nii.gz $1_output_SFNR
fslstats $1_output_CoV.nii.gz -M > $1_output_mean_CoV.txt
fslstats $1_output_SFNR.nii.gz -M > $1_output_mean_SFNR.txt

