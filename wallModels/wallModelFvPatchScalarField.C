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

Class
wallModel

Description
    Base abstract class for LES wall models.

Authors
    Timofey Mukha. All rights reserved.

SourceFiles
    wallModel.C

\*---------------------------------------------------------------------------*/

#include "wallModelFvPatchScalarField.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "turbulenceModel.H"


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
            << "Invalid wall model specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void wallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
}

void wallModelFvPatchScalarField::createFields() const
{
   
    // Create and register h field, if not there already
    if (!db().found("h"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    db().time().timeName(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    

    
    // Create and register uTau field, if not there already.
    // This holds uTau as computed by the wall model.
    if (!db().found("uTau"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "uTau",
                    db().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("uTau", dimVelocity, 0.0),
                h.boundaryField().types()
            )
        );
    }

    if (!db().found("samplingCells"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "samplingCells",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("samplingCells", dimless,0),
                h.boundaryField().types()
            )
        );
    }
    

    // For debugging, create a field that stores uTau as computed by
    // The built-in Spalding law wall model.
    if ((!db().found("uTauBench")) && (debug > 1))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "uTauBench",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("uTauBench", dimVelocity, 0.0),
                h.boundaryField().types()
            )
        );
    }

    if (!db().found("magGradU"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "magGradU",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("magGradU", dimVelocity/dimTime, 0.0),
                h.boundaryField().types()
            )
        );
    }
    
}

tmp<scalarField> wallModelFvPatchScalarField::calcUTauBench
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    scalar kappa_ = 0.4;
    scalar E_ = 9.02501349943;
    
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();
  
    
    forAll(uTau, faceI)
    {
        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[faceI]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[faceI]/nuw[faceI]
                    + magUp[faceI]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[faceI]/nuw[faceI]
                  + magUp[faceI]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[faceI] = max(0.0, ut);

        }
    }
    
    if (db().found("uTauBench"))
    {
    volScalarField & uTauField = const_cast<volScalarField &>(
                                    db().lookupObject<volScalarField>("uTauBench"));
     
    uTauField.boundaryField()[patch().index()] == uTau;
    }
    return tuTau;
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
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
    
    checkType();
    createFields();
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
    cellIndexList_(ptf.cellIndexList_),
    h_(ptf.h_)
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField "
            << "from copy, fvPatch, DimensionedField, and fvPatchFieldMapper"
            << " for patch " << patch().name() << nl;
    }

    checkType();
    //createCellIndexList();
    
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    cellIndexList_(patch().size()),
    //h_(patch().size(), dict.lookupOrDefault<scalar>("h", 0))
    h_(patch().size(), 0)    
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    
    checkType();
    createFields();
    createCellIndexList();  
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf
)
:
    fixedValueFvPatchScalarField(wmpsf),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(wmpsf.h_)
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField "
            << "from copy for patch " << patch().name() << nl;           
    }

    checkType();
    //createCellIndexList();
}


wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wmpsf, iF),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(wmpsf.h_)
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    
    checkType();
    //createCellIndexList();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wallModelFvPatchScalarField::createCellIndexList()
{
    if (debug)
    {
        Info<<"Building sample cell index list for patch " << patch().name()
            << nl;   
    }
    
    //const label size = patch().size();
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField> ("h"));
    
    h_ = h.boundaryField()[patchIndex];

    
    //labelList testCellIndexList(size);

    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh();

    // Create a searcher for the mesh
    meshSearch ms(mesh);
    
    // Grab face centres, normal and adjacent cells' centres
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the indices of adjacent cells
    const labelUList & faceCells = patch().faceCells();

    vector point;
    forAll(faceCentres, i)
    {
        // If h is zero, set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly.
        if (h_[i] == 0)
        {
            if (debug > 1)
            {
                Pout<< "h is 0, for face " << i << " on patch " 
                    << patch().name() << ", using distance to Cn." << nl;   
                
            }
            
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            cellIndexList_[i] = faceCells[i];
                    //ms.findNearestCell(cellCentres[i], -1, true);
            //testCellIndexList[i] = ms.findNearestCell(cellCentres[i], -1, true);
            
        }
        else
        {
            // Grab the point h away along the face normal
            point = faceCentres[i] - faceNormals[i]*h_[i];

            // Check that point is inside the (processor) domain
            // Otherwise fall back to adjacent cell's centre.
            if (!ms.isInside(point))
            {
                if (debug)
                {
                    Pout<< "Point " << point << "is outside the domain. "
                        << "Using Cn." << nl;    
                }

                point = cellCentres[i];
            }

            // Find the cell where the point is located
            cellIndexList_[i] = ms.findNearestCell(point, -1, true);
            
            /*testCellIndexList[i] = ms.findNearestCell(cellCentres[i], -1, true);

            if (cellIndexList_[i] == -1)
            {
                FatalErrorIn
                (
                "void Foam::wallModelFvPatchScalarField::createCellIndexList()\n"
                )   << "Failed to find sampling cell for face " << i << "on patch "
                    << patch().name() << ", with face centre " << faceCentres[i]
                    << abort(FatalError);     

            }*/
            
            // Set h to the distance between face centre and located cell's
            // center
            h_[i] = mag(mesh.C()[cellIndexList_[i]] - faceCentres[i]);
        }
    }
    
    // Assign computed h_ to the global h field
    h.boundaryField()[patch().index()] == h_;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField> ("samplingCells"));
    
    forAll(cellIndexList_, i)
    {
        samplingCells[cellIndexList_[i]] = patchIndex; 
    }
}

void wallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Compute nut and assign
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
