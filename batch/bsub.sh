#!/bin/bash
# https://twiki.cern.ch/twiki/bin/view/Main/BatchJobs#JobSub

jq=8nm
#jq=8nh
#jq=1nd
#jq=2nd

debug=echo
for i in {000..100}; do
  #bsub -R "pool>10000" -q $jq -J $i run.sh $i 10
  bsub -q $jq -J $i run.sh $i 100000
  #sleep 0.1 
  sleep 1
done
