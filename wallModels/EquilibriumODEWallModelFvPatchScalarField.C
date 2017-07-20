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
EquilibriumODEWallModel

Description
    ODE wall model with a zero source term.

Authors
    Timofey Mukha.  All rights reserved.

 * 
\*---------------------------------------------------------------------------*/

#include "EquilibriumODEWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "SpaldingLawOfTheWall.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "dictionary.H"
#include <functional>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void 
EquilibriumODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    ODEWallModelFvPatchScalarField::writeLocalEntries(os);
    //s this is necessary to write entries when decomposing
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


//s this is the constructor used 
EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ODEWallModelFvPatchScalarField(p, iF, dict)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField& wfpsf
)
:
    ODEWallModelFvPatchScalarField(wfpsf)
{}


EquilibriumODEWallModelFvPatchScalarField::
EquilibriumODEWallModelFvPatchScalarField
(
    const EquilibriumODEWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ODEWallModelFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void EquilibriumODEWallModelFvPatchScalarField::write(Ostream& os) const
{
//original:    wallModelFvPatchScalarField::write(os);
//saleh  <<
    ODEWallModelFvPatchScalarField::write(os);
//saleh >>
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    EquilibriumODEWallModelFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
