#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run couette
# Christoph Goniva - Oct. 2014
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
#if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
echo "Building mesh"
blockMesh
decomposePar -force
#fi

# Run CFD
gnome-terminal --title='cfdemSolverPisoIBm couette CFD'  -e "bash parCFDrun.sh" 

