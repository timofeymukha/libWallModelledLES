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

#include "ODEWallModelFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "EddyViscosity.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ODEWallModelFvPatchScalarField, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    eddyViscosity_->write(os);
    os.writeKeyword("eps") << eps_ << token::END_STATEMENT << endl;
    os.writeKeyword("maxIter") << maxIter_ << token::END_STATEMENT << endl;
    os.writeKeyword("nMeshY") << nMeshY_ << token::END_STATEMENT << endl;
}    
    
Foam::scalar 
Foam::ODEWallModelFvPatchScalarField::
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

void Foam::ODEWallModelFvPatchScalarField::createMeshes()
{

    // Number of points in the mesh normal to the wall
     label n=nMeshY_;
           
    forAll(patch(), faceI)
    {
        scalar dx = h_[faceI]/(n -1);

        meshes_[faceI].resize(n, 0.0);
        forAll(meshes_[faceI], pointI)
        {
            // uniform distribution
            meshes_[faceI][pointI] = pointI*dx;
        }
    }
    
}

Foam::tmp<Foam::scalarField>
Foam::ODEWallModelFvPatchScalarField::calcNut() const
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
    
    // Velocity at the boundary
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    
    vectorField Udiff = Uw.patchInternalField() - Uw;
    
    project(Udiff);
    
    // Magnitude of wall-normal velocity gradient
    const vectorField wallGradU(patch().deltaCoeffs()*Udiff);
    
    volVectorField & wallGradUField = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallGradU")
        );
    
    wallGradUField.boundaryField()[patchi] == wallGradU;
    
    // Viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();    
    
    scalarField magGradU = mag(wallGradU);
    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


Foam::tmp<Foam::scalarField>
Foam::ODEWallModelFvPatchScalarField::
calcUTau(const scalarField & magGradU) const
{   
    const label patchi = patch().index();
    const label patchSize = patch().size();
    
    const volScalarField & nuField = db().lookupObject<volScalarField>("nu");
    
    // Velocity and viscosity on boundary
    const fvPatchScalarField & nuw = nuField.boundaryField()[patchi];
    
    // temporary vector for computing the source term
    vector sourceFVec(0, 0, 0);
        
    scalarField magU = mag(U_);
 
    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patchSize, 0.0));
    
    scalarField & uTau = tuTau();            
    
    // Compute uTau for each face
    forAll(uTau, faceI)
    {
        // Points of the 1d wall-normal mesh
        const scalarList & y = meshes_[faceI]; 

        // Starting guess using definition
        scalar tau = (nutw[faceI] + nuw[faceI])*magGradU[faceI];
        
        if (tau > ROOTVSMALL)
        {
            for (int iterI=0; iterI<maxIter_; iterI++)
            {
                scalarList nutValues = 
                    eddyViscosity_->value(y, sqrt(tau), nuw[faceI]);

                scalar integral = integrate(y, 1/(nuw[faceI] + nutValues));
                scalar integral2 = integrate(y, y/(nuw[faceI] + nutValues));
                
                if (mag(integral) < VSMALL )
                {
                    WarningIn
                    (
                        "Foam::ODEWallModelFvPatchScalarField::calcUTau()"
                    )
                        << "when calculating newTau, division by zero occurred."
                        << nl;
                };
                
                
                // Compute the source term
                source(faceI, sourceFVec);
                
                scalar newTau = 
                        sqr(magU[faceI]) + sqr(mag(sourceFVec)*integral2) -
                        2*(U_[faceI] & sourceFVec)*integral2;
                
                newTau  = sqrt(newTau)/integral;
                
                scalar error = mag(tau - newTau)/tau;
                tau = newTau;
                
                if (error < eps_)
                {
                    if (debug > 1)
                    {
                        Info<< "tau_w converged after " << iterI + 1
                            << " iterations." << nl;
                    }
                    break;                            
                }
                

                
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
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("uTauPredicted")
        );
     
    // Assign computed uTau to the boundary field of the global field
    uTauField.boundaryField()[patch().index()] == uTau;
    
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODEWallModelFvPatchScalarField::
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
Foam::ODEWallModelFvPatchScalarField::
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
}


//s this is the constructor when running the code
Foam::ODEWallModelFvPatchScalarField::
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
    nMeshY_(dict.lookupOrDefault<label>("nMeshY", 10))

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
Foam::ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& wfpsf
)
:
    wallModelFvPatchScalarField(wfpsf),
    eddyViscosity_(wfpsf.eddyViscosity_),
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
}


//constructor when running deomposePar/reconstructPar
Foam::ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(wfpsf, iF),
    eddyViscosity_(wfpsf.eddyViscosity_),
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ODEWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}

// ************************************************************************* //
