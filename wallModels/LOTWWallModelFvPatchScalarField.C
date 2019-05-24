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

#include "LOTWWallModelFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "codeRules.H"
#include "scalarListIOList.H"

using namespace std::placeholders;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LOTWWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    rootFinder_->write(os);
    law_->write(os);
}    
    
Foam::tmp<Foam::scalarField> 
Foam::LOTWWallModelFvPatchScalarField::calcNut() const
{
    if (debug)
    {
        Info<< "Updating nut for patch " << patch().name() << nl;        
    }

    const label patchi = patch().index();

    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    
    // Velocity and viscosity on boundary
    const fvPatchScalarField & nuw = nuField.boundaryField()[patchi];

    const scalarListIOList & wallGradU =
        sampler_->db().lookupObject<scalarListIOList>("wallGradU");

    scalarField magGradU(patch().size());
    forAll(magGradU, i)
    {
        magGradU[i] = mag(vector(wallGradU[i][0], wallGradU[i][1], wallGradU[i][2]));
    }

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}

Foam::tmp<Foam::scalarField> 
Foam::LOTWWallModelFvPatchScalarField::
calcUTau(const scalarField & magGradU) const
{  
    const label patchi = patch().index();
    const label patchSize = patch().size();
    
    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    
    // Velocity and viscosity on boundary
    const fvPatchScalarField & nuw = nuField.boundaryField()[patchi];
       
    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patchSize, 0.0));
    scalarField & uTau =
#ifdef FOAM_NEW_TMP_RULES
        tuTau.ref();
#else        
        tuTau();
#endif
    
    // Function to give to the root finder
    std::function<scalar(scalar)> value;
    std::function<scalar(scalar)> derivValue;
    
    // Grab global uTau field
    volScalarField & uTauField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("uTauPredicted")
        );

    // Compute uTau for each face
    forAll(uTau, faceI)
    {
        // Starting guess using old values
        scalar ut = sqrt((nuw[faceI] + nutw[faceI])*magGradU[faceI]);
        
        if (ut > ROOTVSMALL)
        {

            // Construct functions dependant on a single parameter (uTau)
            // from functions given by the law of the wall
            value = std::bind(&LawOfTheWall::value, &law_(), std::ref(sampler_()), faceI,
                              _1, nuw[faceI]);
            
            derivValue = std::bind(&LawOfTheWall::derivative, &law_(),
                                   std::ref(sampler_), faceI, _1, nuw[faceI]);

            // Supply the functions to the root finder
            const_cast<RootFinder &>(rootFinder_()).setFunction(value);
            const_cast<RootFinder &>(rootFinder_()).setDerivative(derivValue);

            // Compute root to get uTau
            uTau[faceI] = max(0.0, rootFinder_->root(ut));

        }
    }
    
    // Assign computed uTau to the boundary field of the global field
#ifdef FOAM_NEW_GEOMFIELD_RULES
    uTauField.boundaryFieldRef()[patchi]
#else        
    uTauField.boundaryField()[patchi]
#endif
    ==
        uTau;
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const fvPatch & p,
    const DimensionedField<scalar, volMesh> & iF
)
:
    wallModelFvPatchScalarField(p, iF),
    rootFinder_(),
    law_(),
    sampler_()
{
    if (debug)
    {
        Info<< "Constructing LOTWwallModelFvPatchScalarField (lotw1) "
            << "from fvPatch and DimensionedField for patch " << patch().name()
            <<  nl;
    }
}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField & ptf,
    const fvPatch & p,
    const DimensionedField<scalar, volMesh> & iF,
    const fvPatchFieldMapper & mapper
)
:
    wallModelFvPatchScalarField(ptf, p, iF, mapper),
    rootFinder_(ptf.rootFinder_, false),
    law_(LawOfTheWall::New(ptf.law_->type(),
                           ptf.law_->constDict())),
    sampler_(new SingleCellSampler(ptf.sampler()))
{
    if (debug)
    {
        Info<< "Constructing LOTWWallModelFvPatchScalarField (lotw2) "
            << "from copy, fvPatch, DimensionedField, and fvPatchFieldMapper"
            << " for patch " << patch().name() << nl;
    }
    law_->addFieldsToSampler(sampler());
}

Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const fvPatch & p,
    const DimensionedField<scalar, volMesh> & iF,
    const dictionary & dict
)
:
    wallModelFvPatchScalarField(p, iF, dict),
    rootFinder_(RootFinder::New(dict.subDict("RootFinder"))),
    law_(LawOfTheWall::New(dict.subDict("Law"))),
    sampler_(new SingleCellSampler(p, averagingTime_))
{
    if (debug)
    {
        Info<< "Constructing LOTWWallModelFvPatchScalarField (lotw3) "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    law_->addFieldsToSampler(sampler());
}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField & wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    rootFinder_(wfpsf.rootFinder_, false),
    law_
    (
        LawOfTheWall::New 
        (
            wfpsf.law_->type(),
            wfpsf.law_->constDict()
        )
    ),
    sampler_(new SingleCellSampler(wfpsf.sampler_()))
{
    if (debug)
    {
        Info<< "Constructing LOTWWallModelFvPatchScalarField (lotw4)"
            << "from copy for patch " << patch().name() << nl;           
    }
}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField & wfpsf,
    const DimensionedField<scalar, volMesh> & iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF),
    rootFinder_(wfpsf.rootFinder_, false),
    law_
    (
        LawOfTheWall::New 
        (
            wfpsf.law_->type(),
            wfpsf.law_->constDict()
        )
    ),
    sampler_(new SingleCellSampler(wfpsf.sampler_()))
{

    if (debug)
    {
        Info<< "Constructing LOTWModelFvPatchScalarField (lotw5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LOTWWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}


void Foam::LOTWWallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    sampler().recomputeFields();
    sampler().sample();

    wallModelFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        LOTWWallModelFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
