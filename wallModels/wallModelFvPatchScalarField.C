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
    defineTypeNameAndDebug(wallModelFvPatchScalarField, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::wallModelFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("wallModelFvPatchScalarField::checkType()")
            << "Invalid wall model specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << nl
            << abort(FatalError);
    }
}


void Foam::wallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("averagingTime")
        << averagingTime_ << token::END_STATEMENT << nl;
}

void Foam::wallModelFvPatchScalarField::createFields() const
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
    /*  if ((!db().found("uTauBench")) && (debug > 1))
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
    */
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    cellIndexList_(patch().size()),
    h_(patch().size(), 0),
    U_(patch().size(), vector(0, 0, 0)),
    averagingTime_(db().time().deltaTValue()) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w1) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
    
    checkType();
    createFields();
    createCellIndexList();
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    cellIndexList_(ptf.cellIndexList_),
    h_(ptf.h_),
    U_(ptf.U_),
    averagingTime_(ptf.averagingTime_) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w2) "
            << "from copy, fvPatch, DimensionedField, and fvPatchFieldMapper"
            << " for patch " << patch().name() << nl;
    }

    checkType();
    //createCellIndexList();   
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    cellIndexList_(patch().size()),
    h_(patch().size(), 0),
    U_(patch().size(), vector(0, 0, 0)),
    averagingTime_(dict.lookupOrDefault<scalar>("averagingTime", 
                                           db().time().deltaTValue()))    
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w3) "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    
    checkType();
    createFields();
    createCellIndexList();  
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf
)
:
    fixedValueFvPatchScalarField(wmpsf),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(wmpsf.h_),
    U_(wmpsf.U_),
    averagingTime_(wmpsf.averagingTime_) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w4)"
            << "from copy for patch " << patch().name() << nl;           
    }

    checkType();
    //createCellIndexList();
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wmpsf, iF),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(wmpsf.h_),
    U_(wmpsf.U_),
    averagingTime_(wmpsf.averagingTime_)     
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    
    checkType();
    //createCellIndexList();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallModelFvPatchScalarField::createCellIndexList()
{
    if (debug)
    {
        Info<<"Building sample cell index list for patch " << patch().name()
            << nl;   
    }
    
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField> ("h"));
    
    h_ = h.boundaryField()[patchIndex];

    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh();

    // Create a searcher for the mesh
    meshSearch ms(mesh);
    
    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the global indices of adjacent cells 
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
                        
            // Set h to the distance between face centre and located cell's
            // center
            h_[i] = mag(mesh.C()[cellIndexList_[i]] - faceCentres[i]);
        }
    }
    
    // Assign computed h_ to the global h field
    h.boundaryField()[patch().index()] == h_;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField> ("samplingCells")
        );
    
    forAll(cellIndexList_, i)
    {
        samplingCells[cellIndexList_[i]] = patchIndex; 
    }
}


void Foam::wallModelFvPatchScalarField::sample()
{
    const volVectorField & UField = db().lookupObject<volVectorField>("U");
    const vectorField & U = UField.internalField();
    const vectorField & UWall = UField.boundaryField()[patch().index()];
    
    // Face normals
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    
    vectorField Up(patch().size());
    
    
    forAll(Up, i)
    {   
        Up[i] = U[cellIndexList_[i]] - UWall[i];
    }
    
    
    scalar eps = db().time().deltaTValue()/averagingTime_;
    
    project(Up);
    
    forAll(U_, i)  
    {    
        U_[i] = eps*Up[i] + (1 - eps)*U_[i];
    }
}

void
Foam::wallModelFvPatchScalarField::project(vectorField & field) const
{
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();

    
    forAll(field, i)
    {    
        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(field[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        field[i] -= normal;        
    }
}

void Foam::wallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    sample();
    
    // Compute nut and assign
    operator==(calcNut());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::wallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //