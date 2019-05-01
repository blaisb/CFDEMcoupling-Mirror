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

#include "DEMbasedDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DEMbasedDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    DEMbasedDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
DEMbasedDrag::DEMbasedDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_))
{
    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implForceDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(7,true); // activate implForceDEMacc switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches(); 

    // check if necessary switches are correct
    if(!forceSubM(0).switches()[2])
        FatalError << "implForceDEM must be set to true!" << abort(FatalError);

    if(!forceSubM(0).switches()[7])
        FatalError << "implForceDEMacc must be set to true!" << abort(FatalError);   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

DEMbasedDrag::~DEMbasedDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DEMbasedDrag::setForce() const
{
    vector position(0,0,0);     // position of the particle
    scalar voidfraction(1);
    vector drag(0,0,0);         // Drag force
    vector dragExplicit(0,0,0);
    vector Ufluid(0,0,0);       // Velocity of the fluid
    label cellI=0;

    int couplingInterval = particleCloud_.dataExchangeM().couplingInterval();

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        drag = vector(0,0,0);
        Ufluid = vector(0,0,0);
        voidfraction = 1.;
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI > -1) // particle Found
        {
            // If on, interpolate the variables to the local position of the particle
            if(forceSubM(0).interpolation())
            {
                position = particleCloud_.position(index);
                voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                Ufluid = UInterpolator_.interpolate(position,cellI);
            }
            //Velocity and voidfraction are assumed to be that of the fluid cell
            else
            {
                voidfraction = voidfraction_[cellI];
                Ufluid = U_[cellI];
            }

            // Drag calculation is carried out within the DEM routine,
            // the CFD one only recuperates the drag vector and sends the volume fraction
            for (int j=0 ; j<3 ; j++) drag[j] = particleCloud_.fAccs()[index][j]/couplingInterval;

            if(forceSubM(0).verbose() && index >=0 && index <2)
            {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Ufluid = " << Ufluid << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "drag = " << drag << endl;
            }

        }        
        // write particle based data to global array using subModel
        forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,1.-voidfraction);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
