# Unit test for the SRFcfdemSolverPiso solver
# This consists in rotating a single paticle in a taylor-Couette
# flow and tracking velocity in the Eulerian Frame of reference
# The profile should be a sinusoidal profile and can be compared with analytical
# solution

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh


#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
results=$casePath/CFD/results

echo "mesh needs to be built"
cd $casePath/CFD
blockMesh > "log_BlockMesh.log"

# Check if the template directory exist, if not create it
if [ -d "$results" ]; then
    echo "Results directory exists, erasing its content"
    rm -r $results
    mkdir -p $results
else
    mkdir $results
fi

echo "Starting simulation with verbosity, this can take up to 5 minutes"


#--------------------------------------------------------------------------------#
#- define variables
logpath=$casePath
headerText="run_parallel_cfdemSolverPiso_couette_CFDDEM"
logDecomposeName="logDecompose_$headerText"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="2"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict          
#--------------------------------------------------------------------------------#

$casePath/parDEMrun.sh > "log_DEM.log"

cd $casePath/CFD
decomposePar -force > "log_Decompose.log"
mpirun -np $nrProcs $solverName -parallel > $casePath/$logfileName
cd ..


echo "Simulation is over, parsing velocity and position data"
# Variables to be parsed from the log output
var1=theta0
var2=thetaf
var3=PhiAngle
var4=radiusParticle
var5=InitialVelocityX
var6=InitialVelocityY
var7=FinalVelocityX
var8=FinalVelocityY
var9=InitialPositionX
var10=InitialPositionY
var11=FinalPositionX
var12=FinalPositionY

# Parse information from the log output using grep, slow but works ok
grep $var1  $casePath/$logfileName | cut -c 18- > $results/out$var1
grep $var2  $casePath/$logfileName | cut -c 18- > $results/out$var2
grep $var3  $casePath/$logfileName | cut -c 18- > $results/out$var3
grep $var4  $casePath/$logfileName | cut -c 18- > $results/out$var4

grep $var5  $casePath/$logfileName | cut -c 18- > $results/out$var5
grep $var6  $casePath/$logfileName | cut -c 18- > $results/out$var6
grep $var7  $casePath/$logfileName | cut -c 18- > $results/out$var7
grep $var8  $casePath/$logfileName | cut -c 18- > $results/out$var8

grep $var9  $casePath/$logfileName | cut -c 18- > $results/out$var9
grep $var10 $casePath/$logfileName | cut -c 18- > $results/out$var10
grep $var11 $casePath/$logfileName | cut -c 18- > $results/out$var11
grep $var12 $casePath/$logfileName | cut -c 18- > $results/out$var12


echo "Parsing is over, plotting results"
cd $casePath/CFD
python3 validateAngle.py


#- clean up case
echo "deleting data at: $casePath ?\n"
read 
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cleanCase
rm $casePath/DEM/post/restart/*
rm $casePath/DEM/post/*
rm $casePath/log*
rm -r $casePath/CFD/results
touch $casePath/DEM/post/.gitignore
touch $casePath/DEM/post/restart/.gitignore
echo "done"
