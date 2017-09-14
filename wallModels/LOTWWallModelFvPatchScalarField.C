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
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "SpaldingLawOfTheWall.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "dictionary.H"
#include <functional>

using namespace std::placeholders;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LOTWWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
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

    // Grab turbulence model to get fields access
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );
    
    // Velocity at the boundary (in case of moving boundary)
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    
    // Magnitude of wall-normal velocity gradient
    const scalarField magGradU(mag(Uw.snGrad()));
    
    // Viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    
    // Debug output
    if ((patch().size() > 0) && (debug))
    {
        // Compute uTau using current model and the default Spalding model
        scalarField uTauNew = calcUTau();
        scalarField uTauBench = calcUTauBench(magGradU);
               
        // Average (not-weighted). As usual, value is local to processor
        scalar avrgNew = sum(uTauNew)/patch().size();
        scalar avrgBench = sum(uTauBench)/patch().size();
        
        // Compute relative error
        scalar diff = mag(avrgNew - avrgBench)/avrgBench*100;
        
        // If error > 1 percent, report
        if (diff > 1)
        {
            Pout<< "Average uTau/uTauBench/diff " << sum(uTauNew)/patch().size() 
                << " " << sum(uTauBench)/patch().size() << " " << diff
                << ", patch " << patch().name() << nl;
        }
    }
    
    return max
    (
        scalar(0),
        sqr(calcUTau())/(magGradU + ROOTVSMALL) - nuw
    );
}

Foam::tmp<Foam::scalarField> 
Foam::LOTWWallModelFvPatchScalarField::calcUTau() const
{

    const label patchi = patch().index();
    const label patchSize = patch().size();
    
    const volVectorField & UField = db().lookupObject<volVectorField>("U");
    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    
    // Velocity in internal field
    const vectorField & U = UField.internalField();

    // Velocity on boundary
    const fvPatchVectorField & Uw = UField.boundaryField()[patchi];
    
    // Magnitude of wall-normal gradient
    const scalarField magGradU(mag(Uw.snGrad()));
    
    volScalarField & gradUField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("magGradU")
        );
    
    gradUField.boundaryField()[patchi] == magGradU;
   
    // Face normals
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    
   
    // Velocity relative to boundary and its magnitude
    vectorField Up(patchSize);
    scalarField magUp(patchSize);
 
    // Normal and parallel components
    vectorField Unormal(patchSize);
    vectorField Upar(patchSize);
    scalarField magUpar(patchSize);
            
    forAll(magUp, i)
    {   
        Up[i] = U[cellIndexList_[i]] - Uw[i];
        magUp[i] = mag(Up[i]);
        
        // Normal component as dot product with (inwards) face normal
        Unormal[i] = -faceNormals[i]*(Up[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        Upar[i] = Up[i] - Unormal[i];
        magUpar[i] = mag(Upar[i]);
    }

    // Viscosity
    const tmp<scalarField> tnuw = nuField.boundaryField()[patchi];
    const scalarField & nuw = tnuw();

    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patchSize, 0.0));
    scalarField& uTau = tuTau();
    
    // Function to give to the root finder
    std::function<scalar(scalar)> value;
    std::function<scalar(scalar)> derivValue;
    
    // Grab global uTau field
    volScalarField & uTauField = 
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("uTau")
        );

    scalarField uTauOld = uTauField.boundaryField()[patchi];

    // Compute uTau for each face
    forAll(uTau, faceI)
    {
        // Starting guess using old values
        scalar ut = sqrt((nuw[faceI] + nutw[faceI])*magGradU[faceI]);

        if (ut > ROOTVSMALL)
        {
            // Construct functions dependant on a single parameter (uTau)
            // from functions given by the law of the wall
            // NOTE we currently still use magUp, not magUpar
            value = std::bind(&LawOfTheWall::value, &law_(), magUp[faceI], 
                              h_[faceI], _1, nuw[faceI]);
            
            derivValue = std::bind(&LawOfTheWall::derivative, &law_(),
                                   magUp[faceI], h_[faceI], _1, nuw[faceI]);

            // Supply the functions to the root finder
            const_cast<RootFinder &>(rootFinder_()).setFunction(value);
            const_cast<RootFinder &>(rootFinder_()).setDerivative(derivValue);
      
            // Compute root to get uTau
            uTau[faceI] = max(0.0, rootFinder_->root(ut));
        }
    }

     
    // Assign computed uTau to the boundary field of the global field
    uTauField.boundaryField()[patchi] == uTau;
    
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(p, iF)
{}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField& ptf,
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
                           ptf.law_->constDict()))
{
}

Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallModelFvPatchScalarField(p, iF, dict),
    rootFinder_(RootFinder::New(dict.subDict("RootFinder"))),
    law_(LawOfTheWall::New(dict.subDict("Law")))
{}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    rootFinder_(wfpsf.rootFinder_),
    law_(wfpsf.law_)
{}


Foam::LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF),
    rootFinder_(wfpsf.rootFinder_),
    law_(wfpsf.law_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LOTWWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
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