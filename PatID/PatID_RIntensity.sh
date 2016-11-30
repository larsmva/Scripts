#!/bin/bash

export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.7/site-packages 

mkdir -p $1/RIntensity
folder="$1/RIntensity"

echo $folder
list=$( find $1 -type f | grep "REGISTERED/GR" | sort)

python $PWD/PatID_RIntensity.py $folder $list
