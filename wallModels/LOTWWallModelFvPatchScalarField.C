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
    Base abstract class for LES wall models.

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

using namespace std::placeholders;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar dummyf(Foam::scalar x){return 0;}

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> LOTWWallModelFvPatchScalarField::calcNut() const
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
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    
    scalarField uTauNew = calcUTau(magGradU);
    scalarField uTauBench = calcUTauBench(magGradU);

    Info<< uTauNew[1] << " " << uTauBench[1] << endl;
      
    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}

tmp<scalarField> LOTWWallModelFvPatchScalarField::calcUTauBench
(
    const scalarField& magGradU
) const
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
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    scalar kappa_ = 0.4;
    scalar E_ = 9.02501349943;
    
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();
  
    
    forAll(uTau, faceI)
    {
        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);
        //Info<< "Face " << faceI << " " << "Inital utau " << ut << endl;
        //Info<< "nut " << nutw[faceI] << " gradU " << magGradU[faceI] << endl; 
        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[faceI]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[faceI]/nuw[faceI]
                    + magUp[faceI]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));
             //   Info<< "Value " << f<<endl;

                scalar df =
                    y[faceI]/nuw[faceI]
                  + magUp[faceI]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;
            //    Info<< "Derivative " << df<<endl;
                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[faceI] = max(0.0, ut);
         //   Info<<uTau[faceI]<< endl;
        }
    }
    return tuTau;
}

tmp<scalarField> LOTWWallModelFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
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
    const scalarField& y = turbModel.y()[patchi];
    
    const vectorField & U = turbModel.U().internalField();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp2(mag(Uw.patchInternalField() - Uw));
    
    scalarField magUp(patch().size());
    
    forAll(magUp, i)
    {
        magUp[i] = mag(U[cellIndexList_[i]] - Uw[i]);
    }
    
   // Info << magUp << nl;
   // Info << magUp2 << nl;
   // Info << "y " << y << nl;
   // Info << "h " << h_ << nl;

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();
       
    std::function<scalar(scalar)> value;
    std::function<scalar(scalar)> derivValue;
        
    dictionary rootFinderDict = dict_.subDict("RootFinder");
    autoPtr<RootFinder> rootFinder = RootFinder::New(dummyf, dummyf, rootFinderDict);
    
    forAll(uTau, faceI)
    {
        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);

        if (ut > ROOTVSMALL)
        {
            value = std::bind(&LawOfTheWall::value, &law_(), magUp[faceI], h_[faceI], _1, nuw[faceI]);
            derivValue = std::bind(&LawOfTheWall::derivative, &law_(), magUp2[faceI], h_[faceI], _1, nuw[faceI]);

            rootFinder->setFunction(value);
            rootFinder->setDerivative(derivValue);
            uTau[faceI] = max(0.0, rootFinder->root(ut));
        }
    }
        
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
    rootFinder_(ptf.rootFinder_)
{}

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
    wallModelFvPatchScalarField(wfpsf)
{}


LOTWWallModelFvPatchScalarField::
LOTWWallModelFvPatchScalarField
(
    const LOTWWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> LOTWWallModelFvPatchScalarField::yPlus() const
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
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void LOTWWallModelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    wallModelFvPatchScalarField::write(os);
    writeLocalEntries(os);
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
