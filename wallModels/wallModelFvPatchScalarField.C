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
#include "wallFvPatch.H"
#include "codeRules.H"


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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    cpuTimeFraction_(0),
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
    const wallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    cpuTimeFraction_(0),
    averagingTime_(ptf.averagingTime_) 
{
    if (debug)
    {
        Info<< "Constructing wallModelFvPatchScalarField (w2) "
            << "from copy, fvPatch, DimensionedField, and fvPatchFieldMapper"
            << " for patch " << patch().name() << nl;
    }

    checkType();   
    createFields();
}


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    cpuTimeFraction_(0),
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


Foam::wallModelFvPatchScalarField::wallModelFvPatchScalarField
(
    const wallModelFvPatchScalarField& wmpsf
)
:
    fixedValueFvPatchScalarField(wmpsf),  
    cpuTimeFraction_(wmpsf.cpuTimeFraction_),
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
    cpuTimeFraction_(wmpsf.cpuTimeFraction_),
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


    scalar startCPUTime = db().time().elapsedCpuTime();
    
    label pI = patch().index();
    
    // Compute nut and assign
    operator==(calcNut());
    
    // Compute uTau
    volVectorField & wss = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallShearStress")
        );
    
    const volVectorField & wallGradUField =
        db().lookupObject<volVectorField>("wallGradU");
    
    const vectorField & wallGradU = wallGradUField.boundaryField()[pI];
    
    const volScalarField & nu = db().lookupObject<volScalarField>("nu");
    const scalarField & nut = *this;

#ifdef FOAM_NEW_GEOMFIELD_RULES
    wss.boundaryFieldRef()[pI]
#else        
    wss.boundaryField()[pI]
#endif
    ==
        (nut + nu.boundaryField()[pI])*wallGradU;

    cpuTimeFraction_ += (db().time().elapsedCpuTime() - startCPUTime);
    Info<< "Wall modelling time consumption = "
        << label(100*cpuTimeFraction_/(db().time().elapsedCpuTime() + SMALL))
        << "%" << nl;
}


void Foam::wallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);  
}


