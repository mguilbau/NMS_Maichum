#!/bin/bash

IN_DIR="/afs/cern.ch/user/m/mguilbau/cumulantStudy/MultiCumulants/"
OUT_DIR="/afs/cern.ch/work/m/mguilbau/cumulant"

jobID=$1
Nevt=$2
if test -z "$jobID"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
 exit 123;
fi
if test -z "$Nevt"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
 exit 123;
fi

echo "Content of working dir folder: "
echo $IN_DIR
ls $IN_DIR

tdir=`mktemp -d`
cd $tdir
cp -r $IN_DIR/Makefile .
cp -r $IN_DIR/bin .
cp -r $IN_DIR/correlations .
cp -r $IN_DIR/data .
cp -r $IN_DIR/env_script.sh .
cp -r $IN_DIR/MultiCumulants .
cp -r $IN_DIR/output .
cp -r $IN_DIR/ToyMC .
cp -r $IN_DIR/vendor .
cp -r $IN_DIR/tests .

echo "Content of tmp dir folder: "
echo $tdir
ls $tdir

source $IN_DIR/env_script.sh
make clean
make
fname="output_toymc_${Nevt}evts_jobID${jobID}"
./bin/toymc.app --generate --analyze --nevents ${Nevt} --harm 2 --output $fname
mkdir -p ${OUT_DIR}/${jobID}
mv $tdir/output/${fname}.root ${OUT_DIR}/${jobID}/.

cd $IN_DIR
rm -rf $tdir
