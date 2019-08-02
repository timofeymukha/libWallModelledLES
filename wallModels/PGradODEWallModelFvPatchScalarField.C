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
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SampledPGradField.H"
#include "scalarListIOList.H"


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
    vectorField & source
) const
{
    // source term = pressure gradient vector projected on the patch face
    
    const scalarListIOList & pGrad =
        sampler_().db().lookupObject<scalarListIOList>("pGrad");

    forAll(source, i)
    {
        source[i] = vector(pGrad[i][0], pGrad[i][1], pGrad[i][2]);
    }
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
    
    sampler_().addField(new SampledPGradField(patch()));
}


Foam::PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const PGradODEWallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ODEWallModelFvPatchScalarField(orig, p, iF, mapper)
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
    
    sampler_().addField(new SampledPGradField(patch()));
}



Foam::PGradODEWallModelFvPatchScalarField::
PGradODEWallModelFvPatchScalarField
(
    const PGradODEWallModelFvPatchScalarField& orig
)
:
    ODEWallModelFvPatchScalarField(orig)
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
    const PGradODEWallModelFvPatchScalarField& orig,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(orig, iF)
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
    ODEWallModelFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        PGradODEWallModelFvPatchScalarField
    );
} 
#endif

// ************************************************************************* //
