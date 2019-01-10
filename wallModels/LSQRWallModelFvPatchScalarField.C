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

#include "LSQRWallModelFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "codeRules.H"
#include "scalarListListIOList.H"
#include "MultiCellSampler.H"
#include "SpaldingLawOfTheWall.H"

using namespace std::placeholders;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LSQRWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    rootFinder_->write(os);
    law_->write(os);
}    
    
Foam::tmp<Foam::scalarField> 
Foam::LSQRWallModelFvPatchScalarField::calcNut() const
{
    if (debug)
    {
        Info<< "Updating nut for patch " << patch().name() << nl;        
    }

    const label patchi = patch().index();

    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    
    // Velocity and viscosity on boundary
    const fvPatchScalarField & nuw = nuField.boundaryField()[patchi];

    const scalarListListIOList & wallGradU =
        sampler_->db().lookupObject<scalarListListIOList>("wallGradU");

    scalarField magGradU(patch().size());
    forAll(magGradU, i)
    {
        magGradU[i] = 
            mag(vector(wallGradU[i][0][0], wallGradU[i][0][1], wallGradU[i][0][2]));
    }

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}

Foam::tmp<Foam::scalarField> 
Foam::LSQRWallModelFvPatchScalarField::
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

    const scalarListList & h = sampler().h();
    const scalarListListList & U = sampler().db().lookupObject<scalarListListIOList>("U");

    // Compute uTau for each face
    forAll(uTau, faceI)
    {

        // Starting guess using old values
        scalar ut = sqrt((nuw[faceI] + nutw[faceI])*magGradU[faceI]);
        
        if (ut > ROOTVSMALL)
        {


            const label n = h[faceI].size();
            scalarList yStar =  h[faceI]*ut/nuw[faceI]; 
            scalarList uStar(n); 

            scalar sumUStar = 0;
            scalar sumLogYStar = 0;
            scalar sumULogYStar = 0;
            scalar sumLogYStar2 = 0;

            forAll(uStar, i)
            {
                uStar[i] =
                    mag(vector(U[faceI][i][0], U[faceI][i][1], U[faceI][i][2]))/ut;

                sumUStar += uStar[i];
                sumLogYStar += log(yStar[i]);
                sumLogYStar2 += sqr(log(yStar[i]));
                sumULogYStar += uStar[i]*log(yStar[i]);
            }

            const scalar kappaNom = n*sumLogYStar2 - sqr(sumLogYStar);
            const scalar kappaDenom = n*sumULogYStar - sumUStar*sumLogYStar;

            const scalar kappa = kappaNom/(kappaDenom + VSMALL);

            const scalar B = (sumUStar - 1/kappa*sumLogYStar)/n;

            Info<< "y* " << yStar << nl;
            Info<< "u* " << uStar << nl;
            Info<< "kappa " << kappa << " B " << B << nl;

            //SpaldingLawOfTheWall law(kappa, B);

            // Construct functions dependant on a single parameter (uTau)
            // from functions given by the law of the wall
            value = std::bind(&LawOfTheWall::value, &law_(), std::ref(sampler_()), faceI,
                              _1, nuw[faceI]);

            //value = std::bind(&SpaldingLawOfTheWall::value, &law, uStar[0]*ut,
                              //yStar[0]*nuw[faceI]/ut, _1, nuw[faceI]);
            
            //derivValue = std::bind(&SpaldingLawOfTheWall::derivative, &law,
                                   //uStar[0]*ut, yStar[0]*nuw[faceI]/ut, _1,
                                   //nuw[faceI]);

            // Supply the functions to the root finder
            //const_cast<RootFinder &>(rootFinder_()).setFunction(value);
            //const_cast<RootFinder &>(rootFinder_()).setDerivative(derivValue);

            // Compute root to get uTau
            //uTau[faceI] = max(0.0, rootFinder_->root(ut));
            uTau[faceI] = ut;

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

Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(p, iF)
{}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(ptf, p, iF, mapper),
    rootFinder_(RootFinder::New(ptf.rootFinder_->type(),
                                ptf.rootFinder_->f(),
                                ptf.rootFinder_->d(),
                                ptf.rootFinder_->eps(),
                                ptf.rootFinder_->maxIter())),
    law_(LawOfTheWall::New(ptf.law_->type(),
                           ptf.law_->constDict())),
    //sampler_(new MultiCellSampler(p, averagingTime_))
    sampler_(new MultiCellSampler(ptf.sampler()))
{
    law_->addFieldsToSampler(sampler());
}

Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallModelFvPatchScalarField(p, iF, dict),
    rootFinder_(RootFinder::New(dict.subDict("RootFinder"))),
    law_(LawOfTheWall::New(dict.subDict("Law"))),
    sampler_(new MultiCellSampler(p, averagingTime_))
{
    law_->addFieldsToSampler(sampler());
}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    rootFinder_(wfpsf.rootFinder_),
    law_(wfpsf.law_),
    sampler_(wfpsf.sampler_)
{}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF),
    rootFinder_(wfpsf.rootFinder_),
    law_(wfpsf.law_),
    sampler_(wfpsf.sampler_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LSQRWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}

void Foam::LSQRWallModelFvPatchScalarField::updateCoeffs()
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
        LSQRWallModelFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
