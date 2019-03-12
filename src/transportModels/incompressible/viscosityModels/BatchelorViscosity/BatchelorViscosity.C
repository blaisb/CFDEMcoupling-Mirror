/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
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

Implementation of the Batchelor viscosity model
\*---------------------------------------------------------------------------*/

#include "BatchelorViscosity.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BatchelorViscosity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        BatchelorViscosity,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BatchelorViscosity::calcNu() const
{

    const volScalarField& voidfraction = U_.mesh().lookupObject<volScalarField>("voidfraction");
    return nu0_*(1.+2.5*(1.-voidfraction)+7.8*(1.-voidfraction)*(1.-voidfraction));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BatchelorViscosity::BatchelorViscosity
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    BatchelorViscosityCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nu0_(BatchelorViscosityCoeffs_.lookup("nu0")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BatchelorViscosity::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    BatchelorViscosityCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");
    BatchelorViscosityCoeffs_.lookup("nu0") >> nu0_;
    
    return true;
}


// ************************************************************************* //
