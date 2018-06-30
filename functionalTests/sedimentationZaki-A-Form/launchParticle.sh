rm DEM/post/*
cd CFD
rm -r 0.*
mpirun -np 2 cfdemSolverPiso -parallel > ../log
#reconstructPar > ../reconstructLog
cd ..

python3 statsZaki.py output
python3 plotSingleCaseZaki.py output

