/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "IBMeshMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IBMeshMotion, 0);

addToRunTimeSelectionTable
(
    meshMotionModel,
    IBMeshMotion,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IBMeshMotion::IBMeshMotion
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    meshMotionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    surfaceName_(propsDict_.lookup("surfaceName")),
    nSurf_(surfaceName_.size()),
    invariant_(propsDict_.lookup("invariant")),
    enableVolumeMeshCorrection_(propsDict_.lookup("enableVolumeMeshCorrection")),
    surf_(new autoPtr<triSurface>[nSurf_]),
    querySurf_(NULL),
    points_(new autoPtr<pointField>[nSurf_]),
    body_
    (
        IOobject
        (
            "body",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        scalar(0)
    ),
    inside_
    (
        IOobject
        (
            "inside",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        scalar(0)
    ),
    cellVolume_
    (
        IOobject
        (
            "cellVolume",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar("zero",dimVolume,0.0)
    ),
    Uo_
    (
        IOobject
        (
            "Uo",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedVector("0", dimensionSet(0,1,-1,0,0),vector(0.,0.,0.))
    ),
    f_
    (
        IOobject
        (
            "fIB",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedVector("0", dimensionSet(0,1,-2,0,0),vector(0.,0.,0.))
    ),
    axis_(propsDict_.lookup("axis")),
    centerOfRotation_(propsDict_.lookup("centerOfRotation")),
    linearVel_(propsDict_.lookup("linearVelocity")),
    omega_(propsDict_.lookup("omega")),
    volumeMesh_(propsDict_.lookup("volumeMesh")),
    omegaV_(vectorList(nSurf_)),
    T_(tensorList(nSurf_)),
    translation_(vectorList(nSurf_)),
    nrParticleCells_(sm.mesh().nPoints()),
    bodyCells_(new autoPtr<labelList>[nSurf_]),
    haloCells_(new autoPtr<labelList>[nSurf_]),
    nBodyCells_(nSurf_),
    nHaloCells_(nSurf_),
    alpha_(readScalar(propsDict_.lookup("alpha"))),
    postProcessing_(propsDict_.lookup("postProcessing")),
    postProcessingFrequency_(readScalar(propsDict_.lookup("postProcessingFrequency"))),
    timePostProcessed_(particleCloud_.mesh().time().value())
{
    //cellVolume_.internalField().field() = &particleCloud_.mesh().V().field(); 
    cellVolume_.ref() = particleCloud_.mesh().V(); 

    // surfaces attribution 
    for (int i=0;i<nSurf_;i++)
    {
        surf_[i].reset(new triSurface(surfaceName_[i]));
        points_[i].reset(new pointField(surf_[i]->points()));

        // Write information about the STL meshes
        surf_[i]->writeStats(Info);

        Info << "About to construct this" << endl;

        // Attribute the memory for the halo and body lists
        haloCells_[i].reset(new labelList(nrParticleCells_));
        //.setSize(nrParticleCells_);
        bodyCells_[i].reset(new labelList(nrParticleCells_));
    }
   
    defineTransformation();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IBMeshMotion::~IBMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IBMeshMotion::defineTransformation()
{
    scalar dt=particleCloud_.mesh().time().deltaT().value();
    
    // ROTATION
    // define omega in rad.s^-1
    for (int i=0; i<nSurf_; i++)
    {
        omegaV_[i] = omega_[i] * axis_[i] * 2.*M_PI/60.;
    }

    // define rotation tensor
    for (int i =0 ; i<nSurf_; i++)
    {
        vector omegaT(axis_[i]*omega_[i]/60.*2.*M_PI*dt);
        vector n1(1,0,0);
        vector n2(Foam::cos(-omegaT[2]),Foam::sin(omegaT[2]),0);
        T_[i]=rotationTensor(n1, n2);
    }

    // TRANSLATION
    for (int i=0 ; i<nSurf_; i++)
    {
    translation_[i] = linearVel_[i] * dt;
    }
}

void IBMeshMotion::moveMesh() const
{
    bool isFound;


    // Lists for halo and for all body
    List<label> bodyCells(nrParticleCells_);
    List<label> haloCells(nrParticleCells_);

    scalar nVel(0);
    vector Up(0,0,0);
    vector lever(0,0,0);

    // Move the mesh
    Info << "Mesh motion beginning : " << particleCloud_.mesh().time().elapsedCpuTime() << endl;
    // Reset body/inside/Uo field
    forAll(body_,cellI)
    {
        body_[cellI]=0.;
        inside_[cellI]=0.;
        Uo_[cellI]= vector(0.,0.,0.);
    }

    // Loop over the surfaces
    for (int i=0 ; i < nSurf_ ; i++)
    {
        // Search operator
        querySurf_.reset(new triSurfaceSearch(surf_[i]));

        // Get new list of vertices inside the mesh
        boolList isInside = querySurf_().calcInside(particleCloud_.mesh().points());

        // Get list of cells inside the mesh
        boolList cellInside = querySurf_().calcInside(particleCloud_.mesh().C());

        // Dual looping to calculate void fraction and establish accessing list
        nBodyCells_[i]=0;
        nHaloCells_[i]=0;
        forAll(particleCloud_.mesh().C(), cellI)
        {
            // Get the list of vertices linked with a cell
            const labelList& vertices = particleCloud_.mesh().cellPoints()[cellI];
            //Reset body field
            nVel=0.;
            isFound=false;

            // Add the body weight due to vertices
            forAll(vertices, vertI)
            {
                body_[cellI] += isInside[vertices[vertI]];
                if (isInside[vertices[vertI]])
                {
                    vector vertexPosition = particleCloud_.mesh().points()[vertices[vertI]];
                    lever = vertexPosition-centerOfRotation_[i];
                    Uo_[cellI] += linearVel_[i] + (omegaV_[i] ^ lever);
                    nVel++;
                    isFound=true;
                }
            }

            // Add the body weight due to the cell centroid
            if(cellInside[cellI])
            {
                body_[cellI] += vertices.size();
                nVel += vertices.size();
                vector cellPosition = particleCloud_.mesh().C()[cellI];
                lever = cellPosition-centerOfRotation_[i];
                Uo_[cellI] +=  vertices.size() * (linearVel_[i] + (omegaV_[i]^lever));
                isFound=true;
            }

            // Rescale since half weight is vertices and half centroid
            // Could be made more generic?
            if (isFound) body_[cellI] = 0.5* body_[cellI]/vertices.size();

            // Generate body and halo list
            if(isFound) 
            {
                bodyCells_[i]()[nBodyCells_[i]]=cellI;
                nBodyCells_[i] ++;
                if (body_[cellI] < 1.)
                {
                    haloCells_[i]()[nHaloCells_[i]]=cellI;
                    nHaloCells_[i]++;
                }
                // Add centroid velocity ponder velocity by the number of velocities
                Uo_[cellI]/= nVel;
                inside_[cellI]=1.;
            }
        }

        if(enableVolumeMeshCorrection_[i])
        correctVolume(i);
        
        // Move mesh points if the mesh is not considered to be static
        if (!invariant_[i])
        {
            pointField temp(points_[i]);
            forAll(temp,pointI)
            {
                //Remove center of rotation before transformation
                temp[pointI] -= centerOfRotation_[i];
            }
            // Rotate
            temp = transform(T_[i],temp); 

            // Regive center of rotation
            forAll(temp,pointI)
            {
                temp[pointI] += centerOfRotation_[i];
            }

            points_[i].reset(new pointField(temp));
            points_[i].reset(new pointField(points_[i]()  +  translation_[i]));

            //Move the center of rotation of the body also
            centerOfRotation_[i] = centerOfRotation_[i] + translation_[i];

            surf_[i]->movePoints(points_[i]()); 
        }

        Info << "Mesh - " << surfaceName_[i] << " - is over" << endl;

    }
    // Pass information to other processor for inside field since it is used in pressure equation
    inside_.correctBoundaryConditions();

    Info << "Mesh motion over : " << particleCloud_.mesh().time().elapsedCpuTime() << endl;
}


void IBMeshMotion::correctVolume(int i) const
{
    // Volume correction utilities
    scalar totalVol;
    scalar haloVol; 
    // Correct body field to have accurate mesh volume
    //-------------------------------------------------
    //First calculate the volume of the body
    // Might be better to do it in a loop actually instead of purely summing

    totalVol=0;
    for (int j =0 ; j < nBodyCells_[i] ; j ++)
    {
        totalVol += cellVolume_[bodyCells_[i]()[j]]*body_[bodyCells_[i]()[j]];
    }

    // Calculate volume already occupied by halo
    haloVol=0.;
    for (int j = 0 ; j < nHaloCells_[i] ; j++)
    {
        haloVol += cellVolume_[haloCells_[i]()[j]]*body_[haloCells_[i]()[j]];
    }

    // Reduce volume accross processors
    reduce(totalVol, sumOp<scalar>());
    reduce(haloVol, sumOp<scalar>());

    // Calculate correction factor

    scalar volFactor = 1. + ((volumeMesh_[i]-totalVol)/haloVol) ;

    if (mag(volFactor) > 8./7.)
    {
        Info << "Vol factor is too large  - " << volFactor << ", it has been saturated" << endl;
        volFactor = 8./7. * volFactor/mag(volFactor);
    }

    //Scale up/down halo to correspond to actual volume
    for (int j =0 ; j < nHaloCells_[i] ; j ++)
    {
        body_[haloCells_[i]()[j]] *= volFactor;
    }

    // temporary statistics for verification
    // ------------------------------------------------------------
    Info << " ------------------------------------------------------" << endl;
    Info << " Volume information for surface : " << surfaceName_[i] << endl;
    Info << " Total volume is : " << totalVol <<endl;
    Info << " Real volume is : " << volumeMesh_[i] << endl;
    Info << " Halo volume local is : " << haloVol << endl;
    Info << " volFactor is : " << volFactor << endl;
    volScalarField bodyVolume(cellVolume_*body_);
    Info << " New total volume is : " << gSum(bodyVolume) << endl;
    Info << "-------------------------------------------------------" << endl;

    // -------------------------------------------------------------
}

void IBMeshMotion::correctF(volVectorField& U) const
{
    f_ = inside_*f_ +  body_ * alpha_/particleCloud_.mesh().time().deltaT() * (Uo_-U);
    f_.correctBoundaryConditions();    
}

void IBMeshMotion::correctUo(volVectorField& U) const
{
    Uo_ = inside_ * ( ((-body_) + 1.) * U + body_ * Uo_);
}



void IBMeshMotion::postProcessing(volScalarField& rho) const
{


    for (int i =0 ; i< nSurf_ ; i++)
    {
        vector drag(0,0,0);
        vector torque(0.,0.,0.);
        vector lever (0.,0.,0.);
        for(int j=0 ; j<nBodyCells_[i] ; j++)
        {
            // Translational force
            drag += cellVolume_[bodyCells_[i]()[j]]*f_[bodyCells_[i]()[j]] * rho[bodyCells_[i]()[j]];
            
            // Angular force
            lever =  particleCloud_.mesh().C()[bodyCells_[i]()[j]]-centerOfRotation_[i];
            torque += particleCloud_.mesh().V()[bodyCells_[i]()[j]] * ((f_[bodyCells_[i]()[j]])^lever) 
                        * rho[bodyCells_[i]()[j]];
        }
        reduce(drag, sumOp<vector>());
        reduce(torque, sumOp<vector>());
        Info << "---> Drag for " << surfaceName_[i] << " is " << drag << endl;
        Info << "---> Torque is " << surfaceName_[i] << " is " << torque << endl;
    
        if (Pstream::master())
        {
            if (particleCloud_.mesh().time().value() >= timePostProcessed_)
            {
                string fName("logForces_"+string(surfaceName_[i].name()));
                std::ofstream log;
                log.open(fName.c_str(),std::ofstream::app);//,std::fstream::out);
                log << particleCloud_.mesh().time().value();
                log << " " << drag[0] << " " << drag[1] << " " << drag[2];
                log << " " << torque[0] << " " << torque[1] << " " << torque[2];
                log << " " << std::endl;
                log.close();
                timePostProcessed_ += postProcessingFrequency_;
            }
        }
    }
}

tmp<volVectorField> IBMeshMotion::setMotion() const
{
    tmp<volVectorField> tsource
        (
         new volVectorField
         (
          IOobject
          (
           "xxx",
           particleCloud_.mesh().time().timeName(),
           particleCloud_.mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE
          ),
          particleCloud_.mesh(),
          dimensionedVector
          (
           "zero",
           dimensionSet(0, 1, -1, 0, 0),
           vector::zero
          )
         )
        );

    moveMesh();

    return tsource;
}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
