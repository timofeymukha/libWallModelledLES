/*---------------------------------------------------------------------------* \
License
    This file is part of libWallModelledLES.

    libWallModelledLES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libWallModelledLES is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with libWallModelledLES. 
    If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SampledPGradField.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Register appropriate fields in the object registry
void Foam::SampledPGradField::registerFields() const
{
    
    // Grab h to copy bcs from it.
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    if (!db().foundObject<volVectorField>("pGrad"))
    {

        volVectorField
        (
            IOobject
            (
                "pGrad",
                db().time().timeName(),
                db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            db(),
            dimensionedVector
            (
                "pGrad",
                dimPressure/dimLength,
                vector(0, 0, 0)
            ),
            h.boundaryField().types()
        );

        
        IOField<vector>
        (
            IOobject
            (
                "pGrad",
                db().time().timeName(),
                db().subRegistry("wallModelSampling", 0),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            vectorField(9, pTraits<vector>::zero)
        );
    }
}

//- Recompute field
void Foam::SampledPGradField::recompute() const
{
    volVectorField & pGrad = const_cast<volVectorField &>
    (
            db().lookupObject<volVectorField>("pGrad")
    );
    const volScalarField & p = db().lookupObject<volScalarField>("p");
    
    pGrad = fvc::grad(p);
    
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //