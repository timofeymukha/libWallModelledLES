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
#include "turbulenceModel.H"
#include "meshSearch.H"
#include "wallFvPatch.H"
#include "codeRules.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(wallModelFvPatchScalarField, 0);
}
#endif

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
    os.writeKeyword("copyToPatchInternalField")
        << copyToPatchInternalField_ << token::END_STATEMENT << nl;
    os.writeKeyword("silent")
        << silent_ << token::END_STATEMENT << nl;
}

void Foam::wallModelFvPatchScalarField::createFields() const
{
    if (debug)
    {
        Info<< "wallModelFvPatchScalarField creating fields" << nl;
    }
    
    // Name of the h field, default to hSampler
    word hName = "hSampler";

    // Check if hSampler exists
    IOobject hHeader
    (
        "hSampler",
        db().time().timeName(),
        db(),
        IOobject::NO_READ
    );
    
    bool foundhSampler = hHeader.typeHeaderOk<volScalarField>();
    db().checkOut("hSampler");

    if (debug)
    {
        Info<< "wallModelFvPatchScalarField: Found hSampler? " << foundhSampler
            << nl;
    }


    if (!foundhSampler)
    {
        Warning
            << "The hSampler field is not found, will try to find h. "
            << "Please note that h will not work with compressible solvers. "
            << "It is recommended to use hSampler in new cases." << nl; 

        IOobject hHeader
        (
            "h",
            db().time().timeName(),
            db(),
            IOobject::NO_READ
        );
        
        if (hHeader.typeHeaderOk<volScalarField>())
        {
            hName = "h";
            if (debug)
            {
                Info<< "wallModelFvPatchScalarField: Found field h" << nl;
            }
        }
        db().checkOut("h");
    }
    
    if (!db().found(hName))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    hName,
                    db().time().timeName(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
      
    const volScalarField & h = db().lookupObject<volScalarField>(hName);
    
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

    // wall-normal gradient field
    if (!db().foundObject<volVectorField>("wallGradU"))
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
                    pTraits<vector>::zero
                ),
                h.boundaryField().types()
            )
        );  
    }
    if (debug)
    {
        Info<< "wallModelFvPatchScalarField finished creating fields" << nl;
    }
}

Foam::tmp<Foam::volScalarField> Foam::wallModelFvPatchScalarField::nu() const
{
    const turbulenceModel& turbModel =
            db().lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
    return turbModel.nu();
}


Foam::tmp<Foam::scalarField> Foam::wallModelFvPatchScalarField::nu
(
    const label patchi
) const
{
    const turbulenceModel& turbModel =
            db().lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
    return turbModel.nu(patchi);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    consumedTime_(0),
    copyToPatchInternalField_(false),
    silent_(false),
    averagingTime_(0)
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w1) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
    
    checkType();
    createFields();
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(orig, p, iF, mapper),
    consumedTime_(0),
    copyToPatchInternalField_(orig.copyToPatchInternalField_),
    silent_(orig.silent_),
    averagingTime_(orig.averagingTime_)
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
    consumedTime_(0),
    copyToPatchInternalField_
    (
        dict.lookupOrDefault<bool>("copyToPatchInternalField", false)
    ),
    silent_(dict.lookupOrDefault<bool>("silent", false)),
    averagingTime_(dict.lookupOrDefault<scalar>("averagingTime", 0))
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w3) "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    
    checkType();
    createFields();
}


#ifdef FOAM_FVPATCHFIELD_NO_COPY
#else
Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& orig 
)
:
    fixedValueFvPatchScalarField(orig),
    consumedTime_(orig.consumedTime_),
    copyToPatchInternalField_(orig.copyToPatchInternalField_),
    silent_(orig.silent_),
    averagingTime_(orig.averagingTime_)
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w4)"
            << "using the copy constructor" << nl;           
    }

    checkType();
}
#endif


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& orig,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(orig, iF),       
    consumedTime_(orig.consumedTime_),
    copyToPatchInternalField_(orig.copyToPatchInternalField_),
    silent_(orig.silent_),
    averagingTime_(orig.averagingTime_)
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


    scalar startCPUTime = db().time().elapsedClockTime();
    
    label pI = patch().index();
    
    // Compute uTau
    volVectorField & wss = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallShearStress")
        );
    
    const volVectorField & wallGradUField =
        db().lookupObject<volVectorField>("wallGradU");
    
    const vectorField & wallGradU = wallGradUField.boundaryField()[pI];
    
    tmp<scalarField> tnuw = this->nu(pI);
    const scalarField& nuw = tnuw();

    // Compute nut and assign
    scalarField nut(calcNut());

    operator==(nut);


    // Assign to the near-wall cells
    if (copyToPatchInternalField())
    {
        volScalarField & nutField = 
            const_cast<volScalarField &>
            (
                db().lookupObject<volScalarField>("nut")
            );
        const labelUList fC = patch().faceCells();

        forAll (nut, i)
        {
            nutField[fC[i]] = nut[i];
        }
    }

#ifdef FOAM_NEW_GEOMFIELD_RULES
    wss.boundaryFieldRef()[pI]
#else        
    wss.boundaryField()[pI]
#endif
    ==
        (nut + nuw)*wallGradU;

    consumedTime_ += (db().time().elapsedClockTime() - startCPUTime);

    // Take the max consumed time across all procs
    reduce(consumedTime_, maxOp<scalar>());
    
    if (!silent_)
    {
        Info<< "Wall modelling time consumption = " << consumedTime_ 
            << "s "  << 100*consumedTime_/(db().time().elapsedClockTime() + SMALL)
            << "% of total " << nl;
    }
}


void Foam::wallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
#ifdef FOAM_NEW_WRITEENTRY
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);  
#endif
}


