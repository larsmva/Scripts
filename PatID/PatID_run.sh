#!/bin/bash
# Job name:
#SBATCH --job-name=kent
#
# Project:
#SBATCH --account=nn9279k
# Wall clock limit:
#SBATCH --time='30:00:00'
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=20000M
#
# Number of tasks (cores):
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=long
##SBATCH --output=output.$SCRATCH 
source /cluster/bin/jobsetup
module load freesurfer
module load fsl

cleanup "mkdir $HOME/$2"
cleanup "cp -r $SCRATCH/MGZ $HOME/$2"
cleanup "cp -r $SCRATCH/REGISTERED $HOME/$2"
cleanup "cp -r $SCRATCH/LTA $HOME/$2"
cleanup "cp -r $SCRATCH/$2 $HOME/$2"

#export PATH=$FREESURFER_HOME/bin:$PATH
echo $SCRATCH
echo `date`

#INPUT $1=filename.tar.gz $2=subjid $3=-bigventricles

tar -xzvf $1 -C $SCRATCH/

mkdir -p $SCRATCH/MGZ
mkdir -p $SCRATCH/REGISTERED
mkdir -p $SCRATCH/LTA


folder=$(find $SCRATCH -type d -links 2 | grep "DICOM")

for f in ${folder} ; do
    ffile=$(find ${f} -type f | head -1) # Finds first Dicom File
    time=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 8 31 | cut -d '.' -f 1)
    date=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 8 21)
    scanseq=$(${FREESURFER_HOME}/bin/mri_probedicom --i ${ffile} --t 18 20)
    mri_convert ${ffile} $SCRATCH/MGZ/${scanseq}-${date}-${time}.mgz
done

find $SCRATCH/MGZ -type f  -size -10M -delete # Delete corrupted files, i.e. files that uses less than 10 MB memory 

T1=$( find $outpath -type f | grep "MGZ/GR" | sort | head -1)
T2=$( find $outpath -type f | grep "MGZ/SE" | sort | head -1)

list=$( find $SCRATCH -type f | grep "MGZ/GR" | sort)

${FREESURFER_HOME}/bin/mri_robust_template --mov ${list} --average 1 --template $SCRATCH/REGISTERED/template.mgz --satit --inittp 1 --fixtp --noit --maxit 10 --subsample 200 
for i in ${list}; do 
	fname=$(basename $i) 
	${FREESURFER_HOME}/bin/mri_robust_register --mov $i --dst $SCRATCH/REGISTERED/template.mgz -satit --maxit 10 --mapmov $$SCRATCH/REGISTERED/$fname --lta $SCRATCH/LTA/"${fname%.*}".lta
done


if [ -z "$T2" ]; then
recon-all -sd $SCRATCH -subjid $2 -i $T1 -all
else 
recon-all -sd $SCRATCH -subjid $2 -i $T1 -T2 $T2 -all
fi



















