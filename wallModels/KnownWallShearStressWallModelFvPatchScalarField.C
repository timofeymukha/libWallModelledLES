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
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "codeRules.H"
#include "scalarListIOList.H"

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
    tmp<scalarField> deltas = patch().deltaCoeffs(); 
    const labelList & faceCells = patch().faceCells(); 


    // Velocity at the boundary
    const volVectorField & UField = db().lookupObject<volVectorField>("U");
    const fvPatchVectorField& Uw = UField.boundaryField()[patchi];
    
    // Magnitude of wall-normal velocity gradient
    scalarField magGradU(patch().size());
    forAll(magGradU, i)
    {
        magGradU[i] = mag((UField[faceCells[i]] - Uw[i])*deltas()[i]);
    }
    
    // Viscosity
    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    const fvPatchScalarField & nuw = nuField.boundaryField()[patchi];
    
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

}

//constructor when running deomposePar/reconstructPar
Foam::KnownWallShearStressWallModelFvPatchScalarField::
KnownWallShearStressWallModelFvPatchScalarField
(
    const KnownWallShearStressWallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(orig, p, iF, mapper)
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
        db().lookupObject<volScalarField>("uTauPredicted")
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
    const KnownWallShearStressWallModelFvPatchScalarField& orig
)
:
    wallModelFvPatchScalarField(orig)
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
    const KnownWallShearStressWallModelFvPatchScalarField& orig,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(orig, iF)
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        KnownWallShearStressWallModelFvPatchScalarField
    );
}
#endif

// ************************************************************************* //
