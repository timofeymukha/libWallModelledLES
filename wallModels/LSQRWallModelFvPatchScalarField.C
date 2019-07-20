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
#include "OFstream.H"

using namespace std::placeholders;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LSQRWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    rootFinder_->write(os);
    law_->write(os);
}    

void Foam::LSQRWallModelFvPatchScalarField::createFields() const
{
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    if (!db().found("kappa"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    db().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar
                (
                    "kappa",
                    dimless,
                    0
                ),
                h.boundaryField().types()
            )
        );
    }

    if (!db().found("B"))
    {
        db().store
        (
            new volScalarField
            (
                IOobject
                (
                    "B",
                    db().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar
                (
                    "B",
                    dimless,
                    0
                ),
                h.boundaryField().types()
            )
        );
    }

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
        //-nuw
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

    // the kappa and B fields to be computed
    scalarField kappa(patchSize, 0.4);
    scalarField B(patchSize, 5.5);
    
    // Function to give to the root finder
    std::function<scalar(scalar)> value;
    std::function<scalar(scalar)> derivValue;
    
    // Grab global uTau field
    volScalarField & uTauField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("uTauPredicted")
        );

    volScalarField & kappaField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("kappa")
        );

    volScalarField & BField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("B")
        );

    const scalarListList & h = sampler().h();
    const scalarListListList & U = sampler().db().lookupObject<scalarListListIOList>("U");

    fileName outputFile(db().time().path()/db().time().timeName()/"lsqr");
    OFstream os(outputFile);
    os << "# n u y kappa B utau" << nl;
 
    // Compute uTau for each face
    forAll(uTau, faceI)
    {
        // Starting guess using old values
        scalar ut = sqrt((nuw[faceI] + nutw[faceI])*magGradU[faceI]);
        
        ut = 0.0414872;
        const label n = h[faceI].size();
        
        if (ut > ROOTVSMALL)
        {
            for (int iterI=0; iterI<35; iterI++)
            {
                const scalarList yStar =  h[faceI]*ut/nuw[faceI]; 
                scalarList u(n);
                forAll(u, i)
                {
                    u[i] =
                        mag(vector(U[faceI][i][0], U[faceI][i][1], U[faceI][i][2]));
                }

                // Need at least 2 points for a linear fit
                //if (n > 1)
                //{


                // HARDCODED THAT 1st CELL IS UNUSED!
                scalarList yStarNew(n-1);
                scalarList uStarNew(n-1);


                scalarList uStar = u/ut; 
                for (int j=0; j<n-1; j++)
                {
                    yStarNew[j] = yStar[j+1];
                    uStarNew[j] = uStar[j+1];
                }

                    scalar sumUStar = sum(uStarNew);
                    scalar sumLogYStar = sum(log(yStarNew));
                    scalar sumULogYStar = sum(uStarNew*log(yStarNew));
                    scalar sumLogYStar2 = sum(sqr(log(yStarNew)));

                // HARDCODED THAT 1st CELL IS UNUSED!
                    const scalar kappaNom = (n -1)*sumLogYStar2 - sqr(sumLogYStar);
                    const scalar kappaDenom = (n -1)*sumULogYStar - sumUStar*sumLogYStar;

                    kappa[faceI] = kappaNom/(kappaDenom + VSMALL);
                    //kappa[faceI]  = max(0.05, min(10, kappa[faceI]));

                    B[faceI]  = (sumUStar - 1/(kappa[faceI]  + SMALL)*sumLogYStar)/(n-1);

                    //Info<< "y* " << yStar << nl;
                    //Info<< "u* " << uStar << nl;
                //}

                if (faceI == 50)
                {
                    Info<< "ut" << ut << nl;
                    Info<< "kappa " << kappa[faceI]  << " B " << B[faceI]  << nl;
                    //Info<< "n" << n << nl;
                    //Info<< "h" << h[faceI] << nl;
                    //Info<< "yStar" << yStar << nl;
                    //Info<< "uStar" << uStar << nl;
                    //Info<< "yStarNew" << yStarNew << nl;
                    //Info<< "uStarNew" << uStarNew << nl;
                    os << n  << " ";
                    forAll (u, uI)
                    {
                        os << u[uI] << " ";
                    }
                    forAll (u, uI)
                    {
                        os << h[faceI][uI] << " ";
                    }
                    os << kappa[faceI] << " ";
                    os << B[faceI] << " ";
                    os << ut << nl;
                }

                SpaldingLawOfTheWall law(kappa[faceI] , B[faceI]);

                // Construct functions dependant on a single parameter (uTau)
                // from functions given by the law of the wall
                value = std::bind
                (
                    static_cast<scalar(SpaldingLawOfTheWall::*)(scalar, scalar, scalar, scalar) const>(&SpaldingLawOfTheWall::value),
                    &law,
                    //u[n-1],
                    //h[faceI][n-1],
                    u[1],
                    h[faceI][1],
                    _1, 
                    nuw[faceI]
                );

                derivValue = std::bind
                (
                    static_cast<scalar(SpaldingLawOfTheWall::*)(scalar, scalar, scalar, scalar) const>(&SpaldingLawOfTheWall::derivative),
                    &law,
                    //u[n-1],
                    //h[faceI][n-1],
                    u[1],
                    h[faceI][1],
                    _1, 
                    nuw[faceI]
                );

                // Supply the functions to the root finder
                const_cast<RootFinder &>(rootFinder_()).setFunction(value);
                const_cast<RootFinder &>(rootFinder_()).setDerivative(derivValue);

                //Info<< ut << nl;

                // Compute root to get uTau
                scalar eps = 0.1;
                eps = 1;
                ut = eps*max(0.0, rootFinder_->root(ut)) + (1-eps)*ut;
            }

            uTau[faceI] = ut;

        }
    }
    
    // Assign computed uTau to the boundary field of the global field
#ifdef FOAM_NEW_GEOMFIELD_RULES
    uTauField.boundaryFieldRef()[patchi] == uTau;
    kappaField.boundaryFieldRef()[patchi] == kappa;
    BField.boundaryFieldRef()[patchi] == B;
#else        
    uTauField.boundaryField()[patchi] == uTau;
    kappaField.boundaryField()[patchi] == kappa;
    BField.boundaryField()[patchi] == B;
#endif
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
{
    createFields();
}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(orig, p, iF, mapper),
    rootFinder_(RootFinder::New(orig.rootFinder_->type(),
                                orig.rootFinder_->f(),
                                orig.rootFinder_->d(),
                                orig.rootFinder_->eps(),
                                orig.rootFinder_->maxIter())),
    law_(LawOfTheWall::New(orig.law_->type(),
                           orig.law_->constDict())),
    //sampler_(new MultiCellSampler(p, averagingTime_))
    sampler_(new MultiCellSampler(orig.sampler()))
{
    law_->addFieldsToSampler(sampler());
    createFields();
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
    createFields();
}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& orig 
)
:
    wallModelFvPatchScalarField(orig),
    rootFinder_
    (
        RootFinder::New 
        (
            orig.rootFinder_->type(),
            orig.rootFinder_->f(),
            orig.rootFinder_->d(),
            orig.rootFinder_->eps(),
            orig.rootFinder_->maxIter()
        )
    ),
    law_
    (
        LawOfTheWall::New 
        (
            orig.law_->type(),
            orig.law_->constDict()
        )
    ),
    sampler_(new MultiCellSampler(orig.sampler_()))
{}


Foam::LSQRWallModelFvPatchScalarField::
LSQRWallModelFvPatchScalarField
(
    const LSQRWallModelFvPatchScalarField& orig,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(orig, iF),
    rootFinder_
    (
        RootFinder::New 
        (
            orig.rootFinder_->type(),
            orig.rootFinder_->f(),
            orig.rootFinder_->d(),
            orig.rootFinder_->eps(),
            orig.rootFinder_->maxIter()
        )
    ),
    law_
    (
        LawOfTheWall::New 
        (
            orig.law_->type(),
            orig.law_->constDict()
        )
    ),
    sampler_(new MultiCellSampler(orig.sampler_()))
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
