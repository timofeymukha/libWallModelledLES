/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
PGradODEWallModel

Description
    ODE wall model with a source term equal to Pressure Gradient.

Authors
    Timofey Mukha.  All rights reserved.

 * 
\*---------------------------------------------------------------------------*/

#include "PGradODEWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "SpaldingLawOfTheWall.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "dictionary.H"
#include <functional>
#include "fvcGrad.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void 
PGradODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    ODEWallModelFvPatchScalarField::writeLocalEntries(os);
    //s this is necessary to write entries when decomposing the domain
}


void 
PGradODEWallModelFvPatchScalarField::source(const int faceI, const vector patchFaceNormal, vector & sourceFVec) const
{
    //pressure gradient vector at cell center associated to h_
    vector pGrad_h=pressGrad_[cellIndexList_[faceI]];

    // Normal component as dot product with (inwards) face normal
    vector pGrad_hNormal = -patchFaceNormal*(pGrad_h & -patchFaceNormal);

    //pressure source term = pressure gradient vector projected on the patch face
    sourceFVec = pGrad_h-pGrad_hNormal;
}

       

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(p, iF)
{
    if (debug)
    {
        Info<< "Constructing PGradODEWallModelFvPatchScalarField (p1) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}


PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const PGradODEWallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ODEWallModelFvPatchScalarField(ptf, p, iF, mapper)
{
    if (debug)
    {
        Info<< "Constructing PGradODEWallModelFvPatchScalarField (p2) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}


//s this is the constructor used 
PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ODEWallModelFvPatchScalarField(p, iF, dict)
{
    if (debug)
    {
        Info<< "Constructing PGradODEWallModelFvPatchScalarField (p3) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}



PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const PGradODEWallModelFvPatchScalarField& wfpsf
)
:
    ODEWallModelFvPatchScalarField(wfpsf)
{
    if (debug)
    {
        Info<< "Constructing PGradODEWallModelFvPatchScalarField (p4) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}


PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const PGradODEWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(wfpsf, iF)
{
    if (debug)
    {
        Info<< "Constructing PGradODEWallModelFvPatchScalarField (p5) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void PGradODEWallModelFvPatchScalarField::write(Ostream& os) const
{
//saleh    wallModelFvPatchScalarField::write(os);
    ODEWallModelFvPatchScalarField::write(os);
}



void PGradODEWallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //calc & update pressGrad_ 
      //look up pressure values for the whole subdomain the patch() belongs to
    const volScalarField& P = db().lookupObject<volScalarField>("p");

      //calculate pressure Gradient
    pressGrad_=fvc::grad(P);

    wallModelFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    PGradODEWallModelFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
