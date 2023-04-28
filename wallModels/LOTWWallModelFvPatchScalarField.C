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
#include "helpers.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "SingleCellSampler.H"

using namespace std::placeholders;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LOTWWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    rootFinder_->write(os);
    law_->write(os);
    sampler_->write(os);
}

Foam::tmp<Foam::scalarField> 
Foam::LOTWWallModelFvPatchScalarField::calcNut() const
{
    if (debug)
    {
        Info<< "Updating nut for patch " << patch().name() << nl;        
    }

    const label patchi = patch().index();

    tmp<scalarField> nuw = this->nu(patchi);
    
    const scalarListIOList & wallGradU =
        sampler_->db().lookupObject<scalarListIOList>("wallGradU");

    scalarField magGradU(Helpers::mag(wallGradU));

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
    
    tmp<scalarField> tnuw = this->nu(patchi);
    const scalarField& nuw = tnuw();
       
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
    sampler_(nullptr)
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
    const LOTWWallModelFvPatchScalarField & orig,
    const fvPatch & p,
    const DimensionedField<scalar, volMesh> & iF,
    const fvPatchFieldMapper & mapper
)
:
    wallModelFvPatchScalarField(orig, p, iF, mapper),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    rootFinder_(orig.rootFinder_.clone()),
    law_(orig.law_.clone()),
#else
    rootFinder_(orig.rootFinder_, false),
    law_(orig.law_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler()))
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
    sampler_
    (
        new SingleCellSampler
        (
            p,
            averagingTime(),
            dict.lookupOrDefault<word>("interpolationType", "cell"),
            dict.lookupOrDefault<word>("sampler", "Tree"),
            dict.lookupOrDefault<word>("lengthScale", "CubeRootVol"),
            dict.lookupOrDefault<bool>("hIsIndex", false)
        )
    )
{
    if (debug)
    {
        Info<< "Constructing LOTWWallModelFvPatchScalarField (lotw3) "
            << "from fvPatch, DimensionedField, and dictionary for patch "
            << patch().name() << nl;
    }
    law_->addFieldsToSampler(sampler());
}


#ifdef FOAM_FVPATCHFIELD_NO_COPY
#else
Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField & orig
)
:
    wallModelFvPatchScalarField(orig),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    rootFinder_(orig.rootFinder_.clone()),
    law_(orig.law_.clone()),
#else
    rootFinder_(orig.rootFinder_, false),
    law_(orig.law_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler_()))
{
    if (debug)
    {
        Info<< "Constructing LOTWWallModelFvPatchScalarField (lotw4)"
            << "from copy for patch " << patch().name() << nl;           
    }
    law_->addFieldsToSampler(sampler());
}
#endif


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField & orig,
    const DimensionedField<scalar, volMesh> & iF
)
:
    wallModelFvPatchScalarField(orig, iF),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    rootFinder_(orig.rootFinder_.clone()),
    law_(orig.law_.clone()),
#else
    rootFinder_(orig.rootFinder_, false),
    law_(orig.law_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler_()))
{

    if (debug)
    {
        Info<< "Constructing LOTWModelFvPatchScalarField (lotw5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    law_->addFieldsToSampler(sampler());
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

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        LOTWWallModelFvPatchScalarField
    );
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
