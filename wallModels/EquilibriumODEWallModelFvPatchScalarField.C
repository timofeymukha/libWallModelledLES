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

#include "EquilibriumODEWallModelFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void 
Foam::EquilibriumODEWallModelFvPatchScalarField::
writeLocalEntries(Ostream& os) const
{
    // This is necessary to write entries when decomposing
    ODEWallModelFvPatchScalarField::writeLocalEntries(os);   
}       

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(p, iF)
{}


Foam::EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ODEWallModelFvPatchScalarField(orig, p, iF, mapper)
{}


Foam::EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh> & iF,
    const dictionary& dict
)
:
    ODEWallModelFvPatchScalarField(p, iF, dict)
{}


Foam::EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField & orig

)
:
    ODEWallModelFvPatchScalarField(orig)
{}


Foam::EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField & orig,
    const DimensionedField<scalar, volMesh> & iF
)
:
    ODEWallModelFvPatchScalarField(orig, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EquilibriumODEWallModelFvPatchScalarField::write(Ostream& os) const
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
        EquilibriumODEWallModelFvPatchScalarField
    );
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

