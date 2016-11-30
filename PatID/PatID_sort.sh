#!/bin/bash

outpath=$2

mkdir -p "${outpath}/MGZ"
mkdir -p "${outpath}/ABEL" 

folder=$(find $1 -type d -links 2)


for f in ${folder} ; do
    ffile=$(find ${f} -type f | head -1) # Finds first Dicom File

    time=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 8 31 | cut -d '.' -f 1)
    date=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 8 21)
    scanseq=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 18 20)
    #Corrupted files are not discarded. 
    mri_convert ${ffile} ${outpath}/MGZ/${scanseq}-${date}-${time}.mgz
    
done

find $outpath/MGZ -type f  -size -10M -delete

