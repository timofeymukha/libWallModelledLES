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

#include "KnownWallShearStressWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "codeRules.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::KnownWallShearStressWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
}    
    

Foam::tmp<Foam::scalarField>
Foam::KnownWallShearStressWallModelFvPatchScalarField::calcNut() const
{
    if (debug)
    {
        Info<< "Updating nut for patch " << patch().name() << nl;        
    }

    
    const label patchi = patch().index();

    // Grab turbulence model to get fields access
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
#ifdef FOAM_NEW_GEOMFIELD_RULES
            internalField().group()
#else        
            dimensionedInternalField().group()
#endif
        )
    );
    
    // Velocity at the boundary
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    
    // Magnitude of wall-normal velocity gradient
    const scalarField magGradU(mag(Uw.snGrad()));
    
    // Viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField & nuw = tnuw();
    
    const volScalarField& tauWallField = db().lookupObject<volScalarField>("tauWall");
    const fvPatchScalarField & tauWall = tauWallField.boundaryField()[patchi];
    
    return max
    (
        scalar(0),
        tauWall/(magGradU + ROOTVSMALL) - nuw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(p, iF)
{

    if (debug)
    {
        Info<< "Constructing KnownWallShearStressWallModelfvPatchScalarField"
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    if (!db().found("tauWall"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "tauWall",
                    db().time().timeName(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    
    label patchI = patch().index();
    
    const volScalarField & tauWall = db().lookupObject<volScalarField>("tauWall");
    volScalarField & uTauField = 
    const_cast<volScalarField &>
    (
        db().lookupObject<volScalarField>("uTau")
    );
   
#ifdef FOAM_NEW_GEOMFIELD_RULES
    uTauField.boundaryFieldRef()[patchI]
#else        
    uTauField.boundaryField()[patchI]
#endif
    ==
        sqrt(tauWall.boundaryField()[patchI]);

}

//constructor when running deomposePar/reconstructPar
Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const KnownWallShearStressWallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(ptf, p, iF, mapper)
{

    if (debug)
    {
        Info<< "Constructing KnownWallShearStressWallModelfvPatchScalarField"
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
}


//s this is the constructor when running the code
Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallModelFvPatchScalarField(p, iF, dict)
{
    if (debug)
    {
        Info<< "Constructing KnownWallShearStressWallModelfvPatchScalarField"
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    if (!db().found("tauWall"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "tauWall",
                    db().time().timeName(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    
    
    label patchI = patch().index();
    
    const volScalarField& tauWall = db().lookupObject<volScalarField>("tauWall");
    
    volScalarField & uTauField = 
    const_cast<volScalarField &>
    (
        db().lookupObject<volScalarField>("uTau")
    );
   
#ifdef FOAM_NEW_GEOMFIELD_RULES
    uTauField.boundaryFieldRef()[patchI]
#else        
    uTauField.boundaryField()[patchI]
#endif
    ==
        sqrt(tauWall.boundaryField()[patchI]);
}


//constructor when running deomposePar/reconstructPar
Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const KnownWallShearStressWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf)
{

    if (debug)
    {
        Info<< "Constructing KnownWallShearStressWallModelfvPatchScalarField"
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
}


//constructor when running deomposePar/reconstructPar
Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const KnownWallShearStressWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF)
{

    if (debug)
    {
        Info<< "Constructing KnownWallShearStressWallModelfvPatchScalarField"
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::KnownWallShearStressWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        KnownWallShearStressWallModelFvPatchScalarField
    );
}


// ************************************************************************* //
