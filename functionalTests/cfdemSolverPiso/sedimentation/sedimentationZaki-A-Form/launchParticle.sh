#Initialize DEM
cd DEM
lmp_auto < in.liggghts_init
rm DEM/post/*
cd ..

#Launch CFD-DEM
cd CFD
rm -r 0.*
blockMesh
decomposePar -force
mpirun -np 2 cfdemSolverPiso -parallel > ../log
#reconstructPar > ../reconstructLog
cd ..

#Post-Process
python3 statsZaki.py output
python3 plotSingleCaseZaki.py output

