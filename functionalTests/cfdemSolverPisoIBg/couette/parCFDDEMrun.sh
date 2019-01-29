#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run cfdemSolverPisoIBm/couette CFD part
# Christoph Goniva - Feb. 2011 ----- Bruno Blais - October 2014
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoIBm_couette_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPisoIBg"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="off"
cleanCase="true"
postproc="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode "true"

# Not sure how you want to integrate this part
python3 validateAngle.py particles.txt
rm particles.txt
python3 validateTorque.py logForces_couette.stl


if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finished? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK
    #source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read
fi

#- clean up case
#if [ $cleanCase == "true" ]
#  then
    echo "deleting data at: $casePath ?\n"
    read
    source $WM_PROJECT_DIR/bin/tools/CleanFunctions
    cd $casePath/CFD
    cleanCase
    rm log*
    cd $casePath
    rm -r $casePath/CFD/clockData
    rm -r $casePath/CFD/particle.txt
    rm -r $casePath/DEM/post/*
    rm -r $casePath/DEM/liggghts.restartCFDEM*
    echo "done"
#fi

#- preserve post directory
touch $casePath/DEM/post/.gitignore

