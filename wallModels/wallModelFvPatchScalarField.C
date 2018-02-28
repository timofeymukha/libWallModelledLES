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
      
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    // Create and register wallShearStress field, if not there already.
    if (!db().found("wallShearStress"))
    {
        db().store
        (
            new volVectorField
            (
                IOobject
                (
                    "wallShearStress",
                    db().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedVector
                (
                    "wallShearStress",
                    sqr(dimVelocity),
                    vector(0, 0, 0)
                ),
                h.boundaryField().types()
            )
        );
    }
    
    // Field with uTau as predicted by the wall model.
    if ((!db().found("uTauPredicted")))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "uTauPredicted",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("uTauPredicted", dimVelocity, 0.0),
                h.boundaryField().types()
            )
        );
    }

    // Wall-normal velocity gradient
    if (!db().found("wallGradU"))
    {
        db().store
        (
            new volVectorField
            (
                IOobject
                (
                    "wallGradU",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedVector
                (
                    "wallGradU",
                    dimVelocity/dimLength,
                    vector(0, 0, 0)
                ),
                h.boundaryField().types()
            )
        );
    }
    
}


void Foam::wallModelFvPatchScalarField::sample()
{
    const volVectorField & Ufield = db().lookupObject<volVectorField>("U");
    const vectorField & U = Ufield.internalField();
    const fvPatchVectorField & Uwall = Ufield.boundaryField()[patch().index()];
      
    // Wall-normal velocity gradient
    vectorField Udiff = Uwall.patchInternalField() - Uwall;
    project(Udiff);
    const vectorField wallGradU(patch().deltaCoeffs()*Udiff);
    
    volVectorField & wallGradUField = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallGradU")
        );
    wallGradUField.boundaryField()[patch().index()] == wallGradU;
    
    // Sampled velocity
    vectorField Up(patch().size()); 
    
    forAll(Up, i)
    {   
        Up[i] = U[cellIndexList_[i]] - Uwall[i];
    }
    
    project(Up);
        
    scalar eps = db().time().deltaTValue()/averagingTime_;
    
    forAll(U_, i)  
    {    
        U_[i] = eps*Up[i] + (1 - eps)*U_[i];
        wallGradU_[i] = eps*wallGradU[i] + (1 - eps)*wallGradU_[i];
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    cellIndexList_(patch()),
    h_(cellIndexList_.h()),
    U_(patch().size(), vector(0, 0, 0)),
    wallGradU_(patch().size(), vector(0, 0, 0)),
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
    
/*    const volVectorField & Ufield = db().lookupObject<volVectorField>("U");
    const vectorField & U = Ufield.internalField();
    const fvPatchVectorField & Uwall = Ufield.boundaryField()[patch().index()];
      
    // Initialize sampled wall-normal velocity gradient
    vectorField Udiff = Uwall.patchInternalField() - Uwall;
    project(Udiff);
    wallGradU_ = patch().deltaCoeffs()*Udiff;
    
    // Initialize sampled velocity velocity gradient
    forAll(U_, i)
    {   
        U_[i] = U[cellIndexList_[i]] - Uwall[i];
    }
*/
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
    h_(cellIndexList_.h()),
    U_(ptf.U_),
    wallGradU_(ptf.wallGradU_),
    averagingTime_(ptf.averagingTime_) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w2) "
            << "from copy, fvPatch, DimensionedField, and fvPatchFieldMapper"
            << " for patch " << patch().name() << nl;
    }

    checkType();   
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    cellIndexList_(patch()),
    h_(cellIndexList_.h()),
    U_(patch().size(), vector(0, 0, 0)),
    wallGradU_(patch().size(), vector(0, 0, 0)),
    averagingTime_
    (
        dict.lookupOrDefault<scalar>
        (
            "averagingTime", db().time().deltaTValue()
        )
    ) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w3) "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    
    checkType();
    createFields();
    
/*    const volVectorField & Ufield = db().lookupObject<volVectorField>("U");
    const vectorField & U = Ufield.internalField();
    const fvPatchVectorField & Uwall = Ufield.boundaryField()[patch().index()];
      
    // Initialize sampled wall-normal velocity gradient
    vectorField Udiff = Uwall.patchInternalField() - Uwall;
    project(Udiff);
    wallGradU_ = patch().deltaCoeffs()*Udiff;
    
    // Initialize sampled velocity velocity gradient
    forAll(U_, i)
    {   
        U_[i] = U[cellIndexList_[i]] - Uwall[i];
    }
*/
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf
)
:
    fixedValueFvPatchScalarField(wmpsf),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(cellIndexList_.h()),
    U_(wmpsf.U_),
    wallGradU_(wmpsf.wallGradU_),        
    averagingTime_(wmpsf.averagingTime_) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w4)"
            << "from copy for patch " << patch().name() << nl;           
    }

    checkType();
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wmpsf, iF),
    cellIndexList_(wmpsf.cellIndexList_),
    h_(cellIndexList_.h()),
    U_(wmpsf.U_),
    wallGradU_(wmpsf.wallGradU_),        
    averagingTime_(wmpsf.averagingTime_)     
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    
    checkType();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::wallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // sample velocity values
    sample();
    
    // Compute nut and assign
    operator==(calcNut());
    
    // Compute uTau
    volVectorField & wss = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallShearStress")
        );
    
    const volVectorField & wallGradU = 
        db().lookupObject<volVectorField>("wallGradU");
    
    const volScalarField & nu = db().lookupObject<volScalarField>("nu");
    
    const scalarField & nut = *this;
    
    label pI = patch().index();
    wss.boundaryField()[pI] == 
        (nut + nu.boundaryField()[pI])*wallGradU.boundaryField()[pI];
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::wallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //