/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"

// Includes for surface mesh operations
#include "searchableSurfaces.H"
#include "pointField.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

// Include to handle multiple meshes
#include "vectorList.H"
#include "tensorList.H"

// Include for dynamic mesh refinement
#include "dynamicFvMesh.H" 
#include "cellSet.H"

#define PI 3.14159265359 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    //#include "createMesh.H"
    //Dynamic mesh refinement
    #include "createDynamicFvMesh.H"
    #include "createTimeControls.H"
    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    #include "declareVars.H"  // declare IBM vars

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Open log file for force and torque
    std::ofstream* dragLog;
    std::ofstream* torqueLog;
    dragLog=new std::ofstream[nSurf];
    torqueLog=new std::ofstream[nSurf];
    for (int i=0 ; i<nSurf ; i++)
    {
        string fName("logDrag_"+string(surfaceName[i].name()));
        dragLog[i].open(fName.c_str());//,std::fstream::out);
        fName="logTorque_"+surfaceName[i].name();
        torqueLog[i].open(fName.c_str());//,std::fstream::out);
    }
    //dragLog.open("logDrag",std::fstream::out);
    //torqueLog.open("logTorque",std::fstream::out);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info << "Dynamic mesh refinement" << endl;

        // Move the mesh
        Info << "Mesh motion beggining : " << runTime.elapsedCpuTime() << endl;
        interFace = mag(mesh.lookupObject<volScalarField>("body"));
        mesh.update(); //dyM
        cellVolume.ref() = mesh.V(); 

        // Reset body/inside/Uo field
        forAll(body,cellI)
        {
            body[cellI]=0.;
            inside[cellI]=0.;
            Uo[cellI]= vector(0.,0.,0.);
        }

        // Loop over the surfaces
        for (int i=0 ; i < nSurf ; i++)
        {
            // Search operator
            querySurf.reset(new triSurfaceSearch(surf[i]));

            // Get new list of vertices inside the mesh
            boolList isInside = querySurf().calcInside(mesh.points());

            // Get list of cells inside the mesh
            boolList cellInside = querySurf().calcInside(mesh.C());


            // Dual looping to calculate void fraction and establish accessing list
            nBodyCells[i]=0;
            nHaloCells[i]=0;
            forAll( mesh.C(), cellI)
            {
                // Get the list of vertices linked with a cell
                const labelList& vertices = mesh.cellPoints()[cellI];
                //Reset body field
                nVel=0.;
                isFound=false;

                // Add the body weight due to vertices
                forAll(vertices, vertI)
                {
                    body[cellI] += isInside[vertices[vertI]];
                    if (isInside[vertices[vertI]])
                    {
                        vector vertexPosition = mesh.points()[vertices[vertI]];
                        lever = vertexPosition-centerOfRotation[i];
                        Uo[cellI] += linearVel[i] + (omegaV[i] ^ lever);
                        nVel++;
                        isFound=true;
                    }
                }

                // Add the body weight due to the cell centroid
                if(cellInside[cellI])
                {
                    body[cellI] += vertices.size();
                    nVel += vertices.size();
                    vector cellPosition = mesh.C()[cellI];
                    lever = cellPosition-centerOfRotation[i];
                    Uo[cellI] +=  vertices.size() * (linearVel[i] + (omegaV[i]^lever));
                    isFound=true;
                }

                // Rescale since half weight is vertices and half centroid
                // Could be made more generic?
                if (isFound) body[cellI] = 0.5* body[cellI]/vertices.size();

                // Generate body and halo list
                if(isFound) 
                {
                    bodyCells[i]()[nBodyCells[i]]=cellI;
                    nBodyCells[i] ++;

                    if (body[cellI] < 1.)
                    {
                        haloCells[i]()[nHaloCells[i]]=cellI;
                        nHaloCells[i]++;
                    }
                    // Add centroid velocity ponder velocity by the number of velocities
                    Uo[cellI]/= nVel;
                    inside[cellI]=1.;
                }
            }


            // Correct body field to have accurate mesh volume
            //-------------------------------------------------
            //First calculate the volume of the body
            // Might be better to do it in a loop actually instead of purely summing
            if (enableVolumeMeshCorrection[i])
            {
                totalVol=0;
                for (int j =0 ; j < nBodyCells[i] ; j ++)
                {
                    totalVol += cellVolume[bodyCells[i]()[j]]*body[bodyCells[i]()[j]];
                }

                // Calculate volume already occupied by halo
                haloVol=0.;
                for (int j = 0 ; j < nHaloCells[i] ; j++)
                {
                    haloVol += cellVolume[haloCells[i]()[j]]*body[haloCells[i]()[j]];
                }

                // Reduce volume accross processors
                reduce(totalVol, sumOp<scalar>());
                reduce(haloVol, sumOp<scalar>());

                // Calculate correction factor

                scalar volFactor = 1. + ((volumeMesh[i]-totalVol)/haloVol) ;

                if (mag(volFactor) > 8./7.)
                {
                    Info << "Vol factor is too large  - " << volFactor << ", it has been saturated" << endl;
                    volFactor = 8./7. * volFactor/mag(volFactor);
                }

                //Scale up/down halo to correspond to actual volume
                for (int j =0 ; j < nHaloCells[i] ; j ++)
                {
                    body[haloCells[i]()[j]] *= volFactor;
                }
            }
            // temporary statistics for verification
            // ------------------------------------------------------------
            /*    Info << " ------------------------------------------------------" << endl;
                Info << " Volume information for surface : " << surfaceName[i] << endl;
                Info << " Total volume is : " << totalVol <<endl;
                Info << " Real volume is : " << volumeMesh[i] << endl;
                Info << " Halo volume local is : " << haloVol << endl;
                Info << " volFactor is : " << volFactor << endl;
                volScalarField bodyVolume(cellVolume*body);
                Info << " New total volume is : " << gSum(bodyVolume) << endl;
                Info << "-------------------------------------------------------" << endl;
            */
            // -------------------------------------------------------------

            // Move mesh points if the mesh is not considered to be static
            if (!invariant[i])
            {
                pointField temp(points[i]);
                //points[i] = transform(T[0], points[0]());
                forAll(temp,pointI)
                {
                    //Remove center of rotation before transformation
                    temp[pointI] -= centerOfRotation[i];
                }
                // Rotate
                temp = transform(T[i],temp); 

                // Regive center of rotation
                forAll(temp,pointI)
                {
                    temp[pointI] += centerOfRotation[i];
                }

                //points[i].reset(new pointField(transform(T[i], points[i]())));
                points[i].reset(new pointField(temp));
                points[i].reset(new pointField(points[i]()  +  translation[i]));
                
                //Move the center of rotation of the body also
                centerOfRotation[i] += translation[i];

                //temp = transform(T[0],temp);
                //temp += translation[i];
                //points[i].reset(new pointField(temp));

                surf[i]->movePoints(points[i]()); 
            }

            Info << " Mesh - " << surfaceName[i] << " - is over" << endl;

        }
        Info << "Mesh motion over : " << runTime.elapsedCpuTime() << endl;
        inside.correctBoundaryConditions();


        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence->divDevReff(U)
        );

        solve(UEqn == -fvc::grad(p) + inside* f);

        // --- PISO loop

        while(piso.correct())
        {
            volScalarField rUA = 1.0/UEqn.A();
            rUA.correctBoundaryConditions();
            U = rUA*UEqn.H();
            phi =  (fvc::interpolate(U) & mesh.Sf())
              + fvc::interpolate(rUA)*fvc::ddtCorr(U, phi);
                
            adjustPhi(phi, U, p);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi) + fvc::div(rUA*inside*f) // used to be - rUA*inside*div(f), but should actually be + 
                );

                pEqn.setReference(pRefCell, pRefValue);
                if (piso.finalNonOrthogonalIter())
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                        pEqn.solve();
                        phi -= pEqn.flux();
                    }
                else
                {
                    pEqn.solve();
                }
            }

            #include "continuityErrs.H"

            U += - rUA*fvc::grad(p) + rUA*inside*f;
            U.correctBoundaryConditions();

            // Define body velocity as a weighted interpolation between Uo and U
            Uo = inside * ( ((-body) + 1.) * U + body * Uo);

            // Correct force field
            f = inside*f +  body * alpha/runTime.deltaT() * (Uo-U);
            f.correctBoundaryConditions();
        }
            
        
        turbulence->correct();
        runTime.write();

        #include "calcForceOnBody.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
