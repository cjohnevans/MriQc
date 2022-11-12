#!/bin/bash
#   quick quality control check of DICOM Siemens data using Tarquin
#   reports CRLB% for tNAA and tCr, FWHM(Hz) and SNR
TARQUINPATH='/cubric/software/tarquin/'

if [ $# -ne 1 ]
then
    echo "mrs_data_check DICOMFILE.dcm"
    exit
fi

if [ ! -f $1 ]
then
    echo "File " $1 " doesn't exist"
    exit
fi

${TARQUINPATH}/tarquin --input $1 --output_txt /tmp/tarquin_results.txt > /dev/null 2>&1
if [ -f /tmp/tarquin.out ]
then
    rm /tmp/tarquin.out
fi

echo $1 > /tmp/tarquin.out
grep 'SNR residual' /tmp/tarquin_results.txt >> /tmp/tarquin.out
grep 'Metab FWHM (Hz)'  /tmp/tarquin_results.txt >> /tmp/tarquin.out
crlb_naa=`grep 'TNAA' /tmp/tarquin_results.txt | awk '{print $3}'`
crlb_cr=`grep 'TCr' /tmp/tarquin_results.txt | awk '{print $3}'`
echo 'CRLB% tNAA       :' $crlb_naa >> /tmp/tarquin.out
echo 'CRLB% tCr        :' $crlb_cr >> /tmp/tarquin.out

cat /tmp/tarquin.out
