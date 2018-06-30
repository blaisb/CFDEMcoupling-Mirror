#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run refineWithinSTL
# Christoph Goniva - Jan. 2015
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi


#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath="$casePath"
headerText="run_refineWithinSTL_utility"
logfileName="log_$headerText"
solverName="refineWithinSTL"
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
#--------------------------------------------------------------------------------#

#- call function to run CFD cas
CFDrun $logpath $logfileName $casePath $headerText $solverName

#- clean up case
echo "deleting data at: $casePath ? otherwise ctrl-C\n"
read
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cd $casePath
cleanCase
echo "done"
