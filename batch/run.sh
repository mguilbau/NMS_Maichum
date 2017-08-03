#!/bin/bash

IN_DIR="/afs/cern.ch/user/m/mguilbau/cumulantStudy/MultiCumulants/"
#OUT_DIR="/afs/cern.ch/work/m/mguilbau/cumulant"
OUT_DIR="/eos/cms/store/user/mguilbau/ToyMC"

jobID=$1
Nevt=$2
Njobs=$3
multmin=$4
multmax=$5
vnfluct=$6
if test -z "$jobID"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
 exit 123;
fi
if test -z "$Nevt"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
 exit 123;
fi
if test -z "$Njobs"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
 exit 123;
fi
if test -z "$multmin"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
 exit 123;
fi
if test -z "$multmax"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
 exit 123;
fi
if test -z "$vnfluct"; then
  echo "Usage as: "
  echo "  - 1st argument: JobID [jobnumbers]"
  echo "  - 2nd argument: N events [nevents]"
  echo "  - 3rd argument: N jobs [njobs]"
  echo "  - 4th argument: minimum multiplicity [multmin]"
  echo "  - 5th argument: maximum multiplicity [multmax]"
  echo "  - 6th argument: w/ or w/o flow fluctuations [vnfluct]"
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
for ijob in `seq 1 ${Njobs}`; do
   fname="output_toymc_${Nevt}evts_jobID${jobID}_mult${multmin}_${multmax}_vnfluct_${vnfluct}_${ijob}"
   ./bin/toymc.app --generate --nevents ${Nevt} --harm 2 --output $fname --multmin ${multmin} --multmax ${multmax} --isVnfluct ${vnfluct}
   mkdir -p ${OUT_DIR}/${jobID}
   mv $tdir/output/${fname}.root ${OUT_DIR}/${jobID}/.
done

cd $IN_DIR
rm -rf $tdir
