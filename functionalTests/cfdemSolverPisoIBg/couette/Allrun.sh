#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run couette
# Christoph Goniva - Oct. 2014
#===================================================================#

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "Building mesh"
    cd $casePath/CFD
    blockMesh
    decomposePar -force
fi

if [ -f "$casePath/DEM/post/liggghts.restart" ];  then
    echo "LIGGGHTS init was run before - using existing restart file"
else
    #- run DEM in new terminal
    $casePath/parDEMrun.sh
fi

#- run parallel CFD-DEM in new terminal
#gnome-terminal --title='cfdemSolverPisoIBm couette CFD'  -e "bash $casePath/parCFDDEMrun.sh" 
bash $casePath/parCFDDEMrun.sh

