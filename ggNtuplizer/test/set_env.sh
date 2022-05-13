#!/bin/bash

cmsenv
echo "cmsenv" executed
source /cvmfs/cms.cern.ch/crab3/crab.sh
echo "crab3" sourced
voms-proxy-init --voms cms --valid 192:00
