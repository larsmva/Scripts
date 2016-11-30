#!/bin/bash

outpath=$1

mkdir -p ${outpath}/REGISTERED
mkdir -p ${outpath}/LTA


list=$( find $outpath -type f | grep "MGZ/GR" | sort)


${FREESURFER_HOME}/bin/mri_robust_template --mov ${list} --average 1 --template ${outpath}/REGISTERED/template.mgz --satit --inittp 1 --fixtp --noit --maxit 10 #--subsample 200 

for i in ${list}; do 
	fname=$(basename $i) 
	${FREESURFER_HOME}/bin/mri_robust_register --mov $i --dst ${outpath}/REGISTERED/template.mgz -satit --maxit 10 --mapmov ${outpath}/REGISTERED/$fname --lta ${outpath}/LTA/"${fname%.*}".lta
done
