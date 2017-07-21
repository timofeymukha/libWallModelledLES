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
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void 
EquilibriumODEWallModelFvPatchScalarField::
writeLocalEntries(Ostream& os) const
{
    // This is necessary to write entries when decomposing
    ODEWallModelFvPatchScalarField::writeLocalEntries(os);   
}       

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(p, iF)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ODEWallModelFvPatchScalarField(ptf, p, iF, mapper)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh> & iF,
    const dictionary& dict
)
:
    ODEWallModelFvPatchScalarField(p, iF, dict)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField & wfpsf
)
:
    ODEWallModelFvPatchScalarField(wfpsf)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField & wfpsf,
    const DimensionedField<scalar, volMesh> & iF
)
:
    ODEWallModelFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void EquilibriumODEWallModelFvPatchScalarField::write(Ostream& os) const
{
    ODEWallModelFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
        fvPatchScalarField,
        EquilibriumODEWallModelFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

