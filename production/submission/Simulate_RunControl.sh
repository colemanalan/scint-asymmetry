#!/bin/bash

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo HERE: $HERE

###########################################
####### Settings
###########################################

ENVSHELL=/home/acoleman/work/icecube-software/surface-array/build/env-shell.sh
FLAGS="--gcd dummy"

###########################################
###########################################

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`

INPUTDIR=$1
OUTPUTDIR=$2
ID=$(($3 + $4))

echo INPUTDIR: $INPUTDIR
echo OUTPUTDIR: $OUTPUTDIR
echo ID: $ID

ID=$(printf "%06d" $ID)
INPUTFILE=$INPUTDIR/$ID/DAT${ID}
echo Input file: $INPUTFILE

if [[ ! -f $INPUTFILE ]]; then
  echo The input file does not exist...quitting
  exit 0
fi

OUTPUTFILE=$OUTPUTDIR/${ID}.i3.gz
TEMPOUTPUTFILE=$OUTPUTDIR/temp_${ID}.i3.gz
echo Output file: $OUTPUTFILE
if [[ -f ${OUTPUTFILE} ]]; then
  echo This file already exists, will not make it again
  exit 1
fi

if [[ ! -d $OUTPUTDIR ]]; then
  echo Making output dir
  mkdir -p $OUTPUTDIR
fi

EXE=$HERE/../SimulateShower.py

echo ""
echo ""
echo Running: $ENVSHELL $EXE $FLAGS -o ${TEMPOUTPUTFILE} $INPUTFILE
echo ""

$ENVSHELL $EXE $FLAGS -o $TEMPOUTPUTFILE $INPUTFILE

if [[ -f $TEMPOUTPUTFILE ]]; then
  echo "Moving the temporary file (temp_${ID}.i3.gz) to its final location $OUTPUTFILE"
  mv $TEMPOUTPUTFILE $OUTPUTFILE
else
  echo "Warning, it seems that the file was not created: $TEMPOUTPUTFILE"
fi

echo "Done with the script!"
