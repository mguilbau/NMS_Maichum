#!/bin/bash

source /Users/jdb/bnl/vendor/root-5.34.36/root/bin/thisroot.sh
source ./env_script.sh

# Generate a random seed for the job
# # crypto quality random number - should also work in parallel since thread safe
trandom="$(od -vAn -N4 -tu4 < /dev/urandom | tr -d '[:space:]')"
echo  "./bin/toymc.app --seed $(($trandom+0)) --generate --harm 2 --nevents 1000 --ptmin 0.3 --ptmax 3.0 --etamin -2.4 --etamax 2.4 --multmin 10 --multmax 20 --isVnfluct 0 --output test"
./bin/toymc.app --seed $(($trandom+0)) --generate --harm 2 --nevents 1000 --ptmin 0.3 --ptmax 3.0 --etamin -2.4 --etamax 2.4 --multmin 10 --multmax 20 --isVnfluct 0 --output test