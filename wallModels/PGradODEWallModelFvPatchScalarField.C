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

#include "PGradODEWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcGrad.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void 
Foam::PGradODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    ODEWallModelFvPatchScalarField::writeLocalEntries(os);
    //s this is necessary to write entries when decomposing the domain
}


void 
Foam::PGradODEWallModelFvPatchScalarField::
source
(
    const int faceI,
    const vector patchFaceNormal,
    vector & sourceFVec
) const
{
    // pressure gradient vector at cell center associated to h_
    vector pGrad = pressGrad_[cellIndexList_[faceI]];

    // Normal component as dot product with (inwards) face normal
    vector pGradNormal = -patchFaceNormal*(pGrad & - patchFaceNormal);

    // source term = pressure gradient vector projected on the patch face
    sourceFVec = pGrad - pGradNormal;
}

       

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PGradODEWallModelFvPatchScalarField::
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


Foam::PGradODEWallModelFvPatchScalarField::
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


Foam::PGradODEWallModelFvPatchScalarField::
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



Foam::PGradODEWallModelFvPatchScalarField::
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


Foam::PGradODEWallModelFvPatchScalarField::
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

void Foam::PGradODEWallModelFvPatchScalarField::write(Ostream& os) const
{
//saleh    wallModelFvPatchScalarField::write(os);
    ODEWallModelFvPatchScalarField::write(os);
}



void Foam::PGradODEWallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //calculate & update pressGrad_ 
    //look up pressure values for the whole subdomain the patch() belongs to
    const volScalarField & P = db().lookupObject<volScalarField>("p");

      //calculate pressure Gradient
    pressGrad_ = fvc::grad(P);

    wallModelFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        PGradODEWallModelFvPatchScalarField
    );
} 

// ************************************************************************* //
