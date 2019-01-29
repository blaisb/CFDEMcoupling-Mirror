/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPisoIBm

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added,
    and additional functionality for IB resolved mesh motion is added.

Contributions
    Bruno Blais
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "OFversion.H"

#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"


#include "cfdemCloud.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh);
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    cfdemCloud particleCloud(mesh);
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.loop())
    {
        Info<< "\nStarting time loop\n" << endl;
            particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            //Smoothen implicit and explicit momCoupling force
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).impParticleForces());
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).expParticleForces());
        }
    
        Info << "update Ksl and f." << endl;
        Ksl = particleCloud.momCoupleM(0).impMomSource();
        Ksl.correctBoundaryConditions();
        f = particleCloud.momCoupleM(1).expMomSource();
        f.correctBoundaryConditions();

       //Force Checks
       //vector fTotal =-gSum( (mesh.V()) * (f.internalField()));
       //reduce(fTotal, sumOp<vector>());
       //vector fImpTotal = gSum(mesh.V()*Ksl.internalField()*(Us.internalField()-U.internalField()));
       //reduce(fImpTotal, sumOp<vector>());
       //Info << "TotalForceExp: " << fTotal << endl;
       //Info << "TotalForceImp: " << fImpTotal << endl;

        #include "solverDebugInfo.H"

        particleCloud.meshMotionM().setMotion();

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
//                + turbulence->divDevReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  - fvm::Sp(Ksl/rho,U)
                );

                UEqn.relax();
                if (piso.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
                    solve(UEqn == - fvc::grad(p) + Ksl/rho*Us + voidfraction*particleCloud.meshMotionM().inside()* particleCloud.meshMotionM().f() - f/rho);
                else if (piso.momentumPredictor())
                    solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us + voidfraction*particleCloud.meshMotionM().inside()* particleCloud.meshMotionM().f() - f/rho);

                // --- PISO loop

                //for (int corr=0; corr<nCorr; corr++)
                //int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);

                while(piso.correct()) //for (int corr=0; corr<nCorrSoph; corr++)
                {
                    volScalarField rUA = 1.0/UEqn.A();
                    rUA.correctBoundaryConditions();
        
                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);

                    U = rUA*UEqn.H();

                    phi = (fvc::interpolate(U*voidfraction) & mesh.Sf() );
                     //+ fvc::ddtPhiCorr(rUAvoidfraction, U, phi);
                    surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                    surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Non-orthogonal pressure corrector loop
                        while (piso.correctNonOrthogonal())
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + particleCloud.ddtVoidfraction() 
                            + fvc::div(
                            rUA*voidfraction*particleCloud.meshMotionM().inside() * particleCloud.meshMotionM().f() *voidfraction
                            -rUA*voidfraction*f/rho
                            )
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                        if (piso.finalNonOrthogonalIter())
                        {
                        // To correct
                            //phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                        }

                        {
                            phiGes -= pEqn.flux();
                            phi=phiGes;
                        }

                    } // end non-orthogonal corrector loop

                    #include "continuityErrorPhiPU.H"
                    
                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA - rUA*particleCloud.meshMotionM().inside()*particleCloud.meshMotionM().f()*voidfraction + rUA*f/rho;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA - rUA*particleCloud.meshMotionM().inside()*particleCloud.meshMotionM().f()*voidfraction + rUA*f/rho;

                    U.correctBoundaryConditions();

                    // Update object velocity to account for halo
                    //particleCloud.meshMotionM().correctUo(U);

                    particleCloud.meshMotionM().correctF(U);
                } // end piso loop
            }

            turbulence->correct();
            particleCloud.meshMotionM().postProcessing(rho);
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
