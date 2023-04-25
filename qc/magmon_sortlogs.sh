#!/bin/bash

# Connectom
dir=`date "+%Y%m%d-%H%M"`
mkdir /cubric/scanners/mri/3t-micro/logs/${dir}
mv /cubric/scanners/mri/3t-micro/logs/*utr /cubric/scanners/mri/3t-micro/logs/${dir}

