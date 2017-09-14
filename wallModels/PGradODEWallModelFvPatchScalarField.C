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
    vector & sourceFVec
) const
{
    // source term = pressure gradient vector projected on the patch face
    sourceFVec = pressureGrad_[faceI];
}


void Foam::PGradODEWallModelFvPatchScalarField::sample()
{
    //calculate & update pressureGrad_ 
    //look up pressure values for the whole subdomain the patch() belongs to
    const volScalarField & p = db().lookupObject<volScalarField>("p");

    //calculate pressure Gradient
    vectorField gradP = fvc::grad(p);
    
    // Face normals
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    
    // pressure gradient vector at cell center associated to h_
    forAll (pressureGrad_, i)
    {
        pressureGrad_[i] = gradP[cellIndexList_[i]];
        vector pGradNormal = -faceNormals[i]*(pressureGrad_[i] & -faceNormals[i]);
        pressureGrad_[i] -= pGradNormal;
    }

    wallModelFvPatchScalarField::sample();
}
       

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(p, iF),
    pressureGrad_(patch().size(), vector(0, 0, 0))
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
    ODEWallModelFvPatchScalarField(ptf, p, iF, mapper),
    pressureGrad_(ptf.pressureGrad_)
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
    ODEWallModelFvPatchScalarField(p, iF, dict),
    pressureGrad_(patch().size(), vector(0, 0, 0))
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
    ODEWallModelFvPatchScalarField(wfpsf),
    pressureGrad_(wfpsf.pressureGrad_)
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
    ODEWallModelFvPatchScalarField(wfpsf, iF),
    pressureGrad_(wfpsf.pressureGrad_)
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
