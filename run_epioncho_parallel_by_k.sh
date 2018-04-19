#!/bin/bash
echo "Assuming R is set up and able to build Rcpp source files"

#  install the parallel package
echo "Checking for GNU parallel"
which parallel
if [[ $? -eq 1 ]]; then
  echo "No GNU parallele. Installing via homebrew"
  brew install parallel
else
  echo "GNU parallel is ready to go"
fi

#  number of cores in the system -1
cores=$(( $(sysctl -n hw.ncpu) - 1 ))
echo "The epioncho R script will run on $cores cores"

#  Get directory of this script no matter where it is called from
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $DIR

#  This file contains the parameters to use
PARAMS="${DIR}/Imperial_parameters.txt"

#  Run Rscripts in parallel
echo "Calling R script. Do not close this terminal window!"
parallel -a ${PARAMS} -j $cores --lb --progress "Rscript $DIR/run_epioncho_ivm_by_k.R {} $DIR {1} {2} {#}"
