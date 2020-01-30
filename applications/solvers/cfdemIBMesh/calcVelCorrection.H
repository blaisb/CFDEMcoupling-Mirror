// Lever arm between center of IB and cell
vector lever(0,0,0);

// calc vel correction
vector Up(0,0,0);
Info << "vel="<< vel << endl;

/*for(int i=0;i<nrParticleCells;i++)
{
    label cellI = particleCells[i];
    lever = U.mesh().C()[cellI]-pos;
    //Info << "cell = " << cellI << endl;
    //scalar rho=1;
    //vector Up = Up + rho*(vel-U[cellI])*mesh.V()[cellI]/mP;
    //U[cellI] = Up+vel;
    U[cellI] = vel + (omega ^ lever) ;
}
U.correctBoundaryConditions();
*/

// Solve pseudo-pressure correction step
fvScalarMatrix phiIBEqn
(
    fvm::laplacian(phiIB)==body*fvc::div(U)
);
solve(phiIBEqn);


//Info << "phiIB = " << phiIB << endl;
//U=U-fvc::grad(phiIB);
//U.correctBoundaryConditions();

// correct the pressure as well 
//p=p+phiIB/runTime.deltaT();  // do we have to  account for rho here? no because P is P/rho in this code
//p.correctBoundaryConditions();


