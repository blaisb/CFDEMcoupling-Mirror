#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
folder=0
template=voidfraction.org
destination=voidfraction
force=postProcessing/forces/0/forces.dat
results=results

# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

# Check if the template directory exist, if not create it
if [ -d "$results" ]; then
    echo "Results directory exists, erasing its content"
    rm -r $results
    mkdir $results
  # Control will enter here if $DIRECTORY exists.
else
    mkdir $results
fi


for i in "1.0" "0.9" "0.8" "0.7" "0.6" "0.5"
do
    cat $folder/$template | sed 's|@@@@@|'$i'|' > $folder/$destination

    echo "voidfraction = " $i
    nonNewtonianVoidfractionIcoFoam > log_$i

    mv "$force" ./$results/$i
done

cd $results
python ../plotViscPlate.py "1.0" "0.9" "0.8" "0.7" "0.6" "0.5"
cd ..

#- clean up case
echo "deleting data at: $casePath ?\n"
read
source $WM_PROJECT_DIR/bin/tools/CleanFunctions
cleanCase
rm -r results
rm log_*
echo "done"
