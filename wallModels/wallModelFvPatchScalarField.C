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
wallModel

Description
    Base abstract class for LES wall models.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    wallModel.C

\*---------------------------------------------------------------------------*/

#include "wallModelFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallModelFvPatchScalarField, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void wallModelFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("wallModelFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void wallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    cellIndexList_(patch().size()),
    h_(patch().size(), 0)
{
    //Info << "From patch and field" << nl;
    checkType();
    createCellIndexList();
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    cellIndexList_(patch().size()),
    h_(patch().size(), 0)
{
    //Info << "From patchField, patch, field and mapper" << nl;
    checkType();
    createCellIndexList();
    
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    dict_(dict),
    cellIndexList_(patch().size()),
    h_(patch().size(), dict.lookupOrDefault<scalar>("h", 0))
{
    //Info << "From patch, field and dict" << nl;
    checkType();
    createCellIndexList();
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    dict_(wfpsf.dict_),
    cellIndexList_(patch().size()),
    h_(patch().size(), 0)
{
    //Info << "From patchField" << nl;
    checkType();
    createCellIndexList();
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    cellIndexList_(patch().size()),
    h_(patch().size())
{
    //Info << "From patchField and field" << nl;
    checkType();
    createCellIndexList();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wallModelFvPatchScalarField::createCellIndexList()
{
    Info<<"Building sample cell index list for patch " << patch().name() << nl;
   
    const label size = patch().size();
    
    labelList testCellIndexList(size);
        
    const fvMesh & mesh = patch().boundaryMesh().mesh();
    
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();

    
    //Info << faceNormals << endl;
    //Info << faceCentres << endl;
    //Info << cellCentres << endl;
    
    vector point;
    forAll(faceCentres, i)
    {
        if (h_[i] == 0)
        {
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
        }
      
        point = faceCentres[i] - faceNormals[i]*h_[i];
        //Info << point << nl;
        cellIndexList_[i] = mesh.findCell(point);
        testCellIndexList[i] = mesh.findCell(cellCentres[i]);
        
        if (cellIndexList_[i] == -1)
        {
            FatalErrorIn
            (
            "void Foam::wallModelFvPatchScalarField::createCellIndexList()\n"
            )   << "Failed to find sampling cell for face " << i
                << " with face centre " << faceCentres[i] << abort(FatalError);     
        }
    }
    
  //  Info << cellIndexList_ << nl;
   // Info << testCellIndexList << nl;
    
}

void wallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(calcNut());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void wallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
