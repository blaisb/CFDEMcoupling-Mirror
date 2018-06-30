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

Class
    IBMeshMotion

SourceFiles
    IBMeshMotion.C

Description
    Realizes moving objects represented by an STL surface using IB.

Acknowledgements
    This code was developped in cooperation with Bruno Blais. 

\*---------------------------------------------------------------------------*/
#ifndef IBMeshMotion_H
#define IBMeshMotion_H

#include "meshMotionModel.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "vectorList.H"
#include "tensorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class IBMeshMotion
:
    public meshMotionModel
{
private:

    dictionary propsDict_;

    // Mesh information
    fileNameList surfaceName_;
    mutable int nSurf_; 
    boolList invariant_;
    boolList enableVolumeMeshCorrection_;
    mutable autoPtr<triSurface>* surf_;
    mutable autoPtr<triSurfaceSearch> querySurf_;
    mutable autoPtr<pointField>* points_;

    // Fields for IBM
    mutable volScalarField body_;
    mutable volScalarField inside_;
    mutable  volScalarField cellVolume_;
    mutable scalarField cellVolumes_;

    mutable volVectorField Uo_;
    mutable volVectorField f_;

    // Motion information
    mutable vectorList axis_; 
    mutable vectorList centerOfRotation_;
    mutable vectorList linearVel_;
    mutable scalarList omega_; // in RPM
    
    mutable scalarList volumeMesh_; 
    mutable vectorList omegaV_;
    mutable tensorList T_;
    mutable vectorList translation_;

    mutable int nrParticleCells_;
    mutable List<label> particleCells_;

    // Halo and body list for force post-processing and volume correction
    mutable autoPtr<labelList>* bodyCells_;
    mutable autoPtr<labelList>* haloCells_;

    mutable List<label> nBodyCells_;
    mutable List<label> nHaloCells_;

    scalar alpha_;

    // Postprocessing information
    bool postProcessing_ ;
    scalar postProcessingFrequency_;
    mutable scalar timePostProcessed_;

    void defineTransformation();
    void correctVolume(int) const;
    void moveMesh() const;



public:

    //- Runtime type information
    TypeName("IBMeshMotion");


    // Constructors

        //- Construct from components
        IBMeshMotion
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~IBMeshMotion();


    // Member Functions
        tmp<volVectorField> setMotion() const;

        void correctF(volVectorField&) const;

        void correctUo(volVectorField&) const;


        virtual tmp<volVectorField> f() const {return f_;};

        virtual tmp<volScalarField> body() const {return body_;};

        virtual tmp<volScalarField> inside() const { return inside_;};

        /*volVectorField& f() const {return f_;};

        volScalarField& body() const {return body_;};

        volScalarField& inside() const {return inside_;};*/

        void postProcessing(volScalarField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
