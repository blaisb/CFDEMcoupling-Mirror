/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 2011-2013 OpenFOAM Foundation
                                Copyright (C) 2009-2012 JKU, Linz
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
    refineWithinSTL

Description
    refines the mesh within a serie of STL meshes and initialize what would 
    be the body velocity

This contribution was developped by
    Bruno Blais, URPEI Ecole Polytechnique Montreal (EPM)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

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
    //Dynamic mesh creation
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    // Declare IBM vars
    #include "declareVars.H"  

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    int count=0;
    Info<< "\nStarting refinement/time loop\n" << endl;
    while(runTime.loop())
    {
        interFace = mag(mesh.lookupObject<volScalarField>("body"));
        mesh.update(); //dyM
       
       // Reset body
        forAll(body,cellI)
        {
            body[cellI]=0.;
        }

        // Loop over the surfaces
        for (int i=0 ; i < nSurf ; i++)
        {
            Info << "Mesh - " << surfaceName[i] << " - is over" << endl;

            // Search operator
            querySurf.reset(new triSurfaceSearch(surf[i]));

            // Get new list of vertices inside the mesh
            boolList isInside = querySurf().calcInside(mesh.points());

            // Get list of cells inside the mesh
            boolList cellInside = querySurf().calcInside(mesh.C());

            // Dual looping to calculate void fraction and establish accessing list
            forAll( mesh.C(), cellI)
            {
                // Get the list of vertices linked with a cell
                const labelList& vertices = mesh.cellPoints()[cellI];

                // Add the body weight due to vertices
                forAll(vertices, vertI)
                {
                    body[cellI] += isInside[vertices[vertI]];
                }

                // Add the body weight due to the cell centroid
                if(cellInside[cellI])
                {
                    body[cellI] += vertices.size();
                }

                // Rescale since half weight is vertices and half centroid
                // Could be made more generic?
                // Body is defined as half so generic dynamicMesh dictionnary can be used
                body[cellI] = 0.25* body[cellI]/vertices.size();
            }
            Info << " Mesh - " << surfaceName[i] << " - is over" << endl;
        }
        Info<< "Time = " << runTime.timeName() << nl << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        count++;
        if (count==5) break;

    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
