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
ODEWallModel

Description
    Class for wall models based on an ODE.

Authors
    Timofey Mukha.  All rights reserved.

 * 
\*---------------------------------------------------------------------------*/

#include "ODEWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "SpaldingLawOfTheWall.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "dictionary.H"
#include <functional>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ODEWallModelFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void ODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    eddyViscosity_->write(os);
    os.writeKeyword("eps") << eps_ << token::END_STATEMENT << endl;
    os.writeKeyword("maxIter") << maxIter_ << token::END_STATEMENT << endl;
    os.writeKeyword("nMeshY") << nMeshY_ << token::END_STATEMENT << endl;
}    
    

scalar ODEWallModelFvPatchScalarField::
integrate(const scalarList & y, const scalarList & v) const
{
    
    // trapezoidal rule for now
    scalar integral = 0;
    
    for (int i=0; i<y.size()-1; i++)
    {
        integral += (y[i+1] - y[i])*(v[i+1] + v[i]);
    }
    
    return 0.5*integral;
}

void ODEWallModelFvPatchScalarField::createMeshes()
{

    // Number of points in the mesh normal to the wall
     label n=nMeshY_;
           
    forAll(patch(), faceI)
    {
        scalar dx = h_[faceI]/(n -1);

        meshes_[faceI].resize(n, 0.0);
        forAll(meshes_[faceI], pointI)
        {
            // uniform distribution for now..
            meshes_[faceI][pointI] = pointI*dx;
        }
    }
    
}

tmp<scalarField> ODEWallModelFvPatchScalarField::calcNut() const
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


tmp<scalarField> ODEWallModelFvPatchScalarField::calcUTau() const
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
    
    // Velocity in internal field
    const vectorField & U = turbModel.U().internalField();

    // Velocity on boundary
    const fvPatchVectorField & Uw = turbModel.U().boundaryField()[patchi];
    
    // Magnitude of wall-normal gradient
    const scalarField magGradU(mag(Uw.snGrad()));

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

    // temporary vectors for computing source term
    vector sourceFVec(0,0,0);
    vector patchFaceNormal(0,0,0);
            
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
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField & nuw = tnuw();

    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();
    
        
    
    // Compute uTau for each face
    forAll(uTau, faceI)
    {
        const scalarList & y = meshes_[faceI]; //auxiliary points normal to the patch

        // Starting guess using definition
        scalar tau = (nutw[faceI] + nuw[faceI])*magGradU[faceI];
        
        //Info << "Guess " << uTau[faceI] << nl;
        if (tau > ROOTVSMALL)
        {
            for (int iterI=0; iterI<maxIter_; iterI++)
            {
                scalarList nutValues = 
                    eddyViscosity_->value(y, sqrt(tau), nuw[faceI]);

                scalar integral = integrate(y, 1/(nuw[faceI] + nutValues));
                scalar integral2 = integrate(y, y/(nuw[faceI] + nutValues));

                 // compute source term
                patchFaceNormal=faceNormals[faceI];
                source(faceI,patchFaceNormal,sourceFVec);
                
                scalar newTauT0 = Upar[faceI][0] - sourceFVec[0] * integral2;
                scalar newTauT2 = Upar[faceI][2] - sourceFVec[2] * integral2;
                scalar newTauT  = sqr(newTauT0) + sqr(newTauT2);

//                scalar newTauT1 = sqr(Upar[faceI][0])+sqr(Upar[faceI][2]);
//                scalar newTauT2 = -2.0*integral2*(Upar[faceI][0]*sourceFVec[0]+Upar[faceI][2]*sourceFVec[2]);
//                scalar newTauT3 = sqr(integral2)*(sqr(sourceFVec[0])+sqr(sourceFVec[2]));
//                scalar newTauT = newTauT1+newTauT2+newTauT3;

//                if (newTauT < 0 ) {
//                   WarningIn("Foam::ODEWallModelFvPatchScalarField::calcUTau()")
//                        << "when calculating newTau, sqrt of a negative value occurred. " << nl;
//                };

                if (integral == 0 ) {
                   WarningIn("Foam::ODEWallModelFvPatchScalarField::calcUTau()")
                        << "when calculating newTau, division by zero occurred. " << nl;
                };

                scalar newTau = sqrt(newTauT)/integral;
                
                scalar error = mag(tau - newTau)/tau;
                
                if (error < eps_)
                {
                    if (debug > 1)
                    {
                        Info<< "tau_w converged after " << iterI
                            << " iterations." << nl;
                    }
                    break;                            
                }
                
                tau = newTau;
                
                if ((debug > 1) && (iterI == maxIter_-1))
                {
                    WarningIn("Foam::ODEWallModelFvPatchScalarField::calcUTau()")
                        << "tau_w did not converge to desired tolerance "
                        << eps_ << ". Error value: " << error << nl;
                }    
            }
            
            uTau[faceI] = max(0.0, sqrt(tau));
        }
    }

    // Grab global uTau field
    volScalarField & uTauField = 
        const_cast<volScalarField &>(db().lookupObject<volScalarField>("uTau"));
     
    // Assign computed uTau to the boundary field of the global field
    uTauField.boundaryField()[patch().index()] == uTau;
    
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(p, iF),
    meshes_(patch().size()),
    maxIter_(10),
    eps_(1e-3),
    nMeshY_(5)
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o1) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    createMeshes();
}

//constructor when running deomposePar/reconstructPar
ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(ptf, p, iF, mapper),
    eddyViscosity_(EddyViscosity::New(ptf.eddyViscosity_->type(),
                   ptf.eddyViscosity_->constDict())),
//    meshes_(patch().size()),  
    meshes_(ptf.meshes_),
    maxIter_(ptf.maxIter_),
    eps_(ptf.eps_),
    nMeshY_(ptf.nMeshY_)
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o2) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    //createMeshes();
}


//s this is the constructor when running the code
ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallModelFvPatchScalarField(p, iF, dict),
    eddyViscosity_(EddyViscosity::New(dict.subDict("EddyViscosity"))),
    meshes_(patch().size()),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 10)),
    eps_(dict.lookupOrDefault<scalar>("eps", 1e-3)),
    nMeshY_(dict.lookupOrDefault<label>("nMeshY", 5))

{
    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o3) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    createMeshes();

}

//constructor when running deomposePar/reconstructPar
ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    eddyViscosity_(wfpsf.eddyViscosity_),
//    meshes_(patch().size()),
    meshes_(wfpsf.meshes_),
    maxIter_(wfpsf.maxIter_),
    eps_(wfpsf.eps_),
    nMeshY_(wfpsf.nMeshY_)
    
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o4) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

//    createMeshes();
}


//constructor when running deomposePar/reconstructPar
ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF),
    eddyViscosity_(wfpsf.eddyViscosity_),
//    meshes_(patch().size()),
    meshes_(wfpsf.meshes_),
    maxIter_(wfpsf.maxIter_),
    eps_(wfpsf.eps_),
    nMeshY_(wfpsf.nMeshY_)
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

//    createMeshes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ODEWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}


} // End namespace Foam

// ************************************************************************* //
