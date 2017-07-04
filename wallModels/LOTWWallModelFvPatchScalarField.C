/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
LOTWWallModel

Description
    Class for wall models based on a Law of the Wall.

Authors
    Timofey Mukha.  All rights reserved.

 * 
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

//saleh, for gradient
#include "fvcGrad.H"


using namespace std::placeholders;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void LOTWWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    rootFinder_->write(os);
    law_->write(os);
}    
    
tmp<scalarField> LOTWWallModelFvPatchScalarField::calcNut() const
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



tmp<scalarField> LOTWWallModelFvPatchScalarField::calcUTau() const
{

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );
    //const scalarField& y = turbModel.y()[patchi];
    
    // Velocity in internal field
    const vectorField & U = turbModel.U().internalField();

    // Velocity on boundary
    const fvPatchVectorField & Uw = turbModel.U().boundaryField()[patchi];
    
    // Magnitude of wall-normal gradient
    const scalarField magGradU(mag(Uw.snGrad()));
    volScalarField & gradUField = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField>("magGradU"));
    gradUField.boundaryField()[patch().index()] == magGradU;
   
    // Face normals
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    
   
    // Velocity relative to boundary and its magnitude
    vectorField Up(patch().size());
    scalarField magUp(patch().size());
 
    // Normal and parallel components
    vectorField Unormal(patch().size());
    vectorField Upar(patch().size());
    scalarField magUpar(patch().size());




    //saleh start >>>>>>>>>>>>>>>>>>>>>>>
//    const scalarField & P = turbModel.p().internalField();

    //look up pressure values
    const volScalarField& P = db().lookupObject<volScalarField>("p");
//    const volVectorField& UU = db().lookupObject<volVectorField>("U");


    Pout << "pout size of P = " << P.size() <<nl;
//    Info << "Info size of P = " << P.size() <<nl;
    Pout << "P[0] = " << P[0] <<nl;

      volVectorField GradP = fvc::grad(P);
//    volTensorField GradU = fvc::grad(UU);


//    Pout << "gradP[0]=" << GradP[0] << nl;
//  vectorField& GradP = P.fvc::grad();
//    volTensorField gradU = fvc::grad(U);


    //Saleh end  <<<<<<<<<<<<<<<<<<<<<<<<<

            
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
    
    /*Info << "Normals " << faceNormals << nl << nl;;
    Info << Up << nl << nl;
    Info << "Unormal " << Unormal << nl << nl;
    Info << "Upar " << Upar << nl;
    Info << "diff " << magUp - sqrt((sqr(magUpar) + sqr(mag(Unormal)))) << nl;
    */

    // Viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField & nuw = tnuw();

    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();
    
    // Function to give to the root finder
    std::function<scalar(scalar)> value;
    std::function<scalar(scalar)> derivValue;
    
    // Grab global uTau field
    volScalarField & uTauField = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField>("uTau"));

    scalarField uTauOld = uTauField.boundaryField()[patch().index()];

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
    uTauField.boundaryField()[patch().index()] == uTau;
    
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(p, iF)
{}


LOTWWallModelFvPatchScalarField::
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
{}


//s this is the contsructor
LOTWWallModelFvPatchScalarField::
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


LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    rootFinder_(wfpsf.rootFinder_),
    law_(wfpsf.law_)
{}


LOTWWallModelFvPatchScalarField::
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

void LOTWWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    LOTWWallModelFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
