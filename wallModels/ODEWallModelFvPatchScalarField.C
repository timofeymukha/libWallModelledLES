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
#include "addToRunTimeSelectionTable.H"
#include "codeRules.H"
#include "scalarListIOList.H"
#include "helpers.H"
#include "AdaptiveIntegrator.hpp"
#include "EquilibriumODEExplicitLawOfTheWall.H"
#include <functional>

typedef std::function<Foam::scalar(const Foam::scalar)> IntegrandFunc;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(ODEWallModelFvPatchScalarField, 0);
}
#endif


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ODEWallModelFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    wallModelFvPatchScalarField::writeLocalEntries(os);
    eddyViscosity_->write(os);
    sampler_->write(os);
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

    if (debug)
    {
        Info<< "Creating 1D mesh for patch " << patch().name();
    }

    // Number of points in the mesh normal to the wall
     label n=nMeshY_;

    forAll(patch(), faceI)
    {
        scalar dx = sampler().h()[faceI]/(n -1);

        meshes_[faceI].resize(n, 0.0);
        forAll(meshes_[faceI], pointI)
        {
            // uniform distribution
            meshes_[faceI][pointI] = pointI*dx;
        }
    }

    if (debug)
    {
        Info<< " Done" << nl;;
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
Foam::ODEWallModelFvPatchScalarField::
calcUTau(const scalarField & magGradU) const
{

    const label patchi = patch().index();
    const label patchSize = patch().size();


    tmp<scalarField> tnuw = this->nu(patchi);
    const scalarField& nuw = tnuw();
    
    // Grab iterations field
    auto & iterationsField =
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("solverIterations")
        );

    // Grab wm info
    auto & wallModellingInfoField =
        const_cast<volTensorField &>
        (
            db().lookupObject<volTensorField>("wallModellingInfo")
        );
    
    // Iterations of the non-linear solver
//    tmp<vectorField> tIterations(new vectorField(patchSize, vector(0,0,0)));
//    vectorField & iterations = tIterations.ref();
    vectorField iterations = vectorField(patchSize, vector(0,0,0));
    iterations = vector(0, 0, 0);

    // Store diagonstics stuff about explicit wall modelling
    tmp<tensorField> twallModellingInfo(new  tensorField(patchSize));
    tensorField & wallModellingInfo = twallModellingInfo.ref();
    wallModellingInfo = wallModellingInfoField.boundaryFieldRef()[patchi];

    // vectorField for storing the source term
    vectorField sourceField(patchSize, vector(0, 0, 0));

    // Compute the source term
    source(sourceField);

    const scalarListIOList & U = sampler().db().lookupObject<scalarListIOList>("U");
    scalarField magU(patch().size());

    forAll(magU, i)
    {
        magU[i] = mag(vector(U[i][0], U[i][1], U[i][2]));
    }

    // Turbulent viscosity
    const scalarField & nutw = *this;

    // Computed uTau
    tmp<scalarField> tuTau(new scalarField(patchSize, 0.0));

    scalarField & uTau = tuTau.ref();

    AdaptiveIntegrator<scalar (scalar)> quad;

    scalar kappa = 0.387;
    scalar Aplus = 15.2516;

    EquilibriumODEExplicitLawOfTheWall EqOde(kappa, Aplus);

    if (patchSize)
    {
//        Pout << "Rrank " << Pstream::myProcNo() << " " << max(iterations) << endl;
    }


    // Compute uTau for each face
    forAll(uTau, faceI)
    {

        // Points of the 1d wall-normal mesh
        const scalarList & y = meshes_[faceI];

        // Starting guess using definition
        scalar tau = (nutw[faceI] + nuw[faceI])*magGradU[faceI];

        // Starting guess using old values
        if (wallModellingInfo[faceI][0] > VSMALL)
        {
//            tau = wallModellingInfo[faceI][0];
        }
        scalar old_tau = tau;


        scalar utau_exp = EqOde.uTau(sampler_, faceI, nuw[faceI]);

        wallModellingInfo[faceI] = tensor(0,0,0,0,0,0,0,0,0);

        //if (tau > ROOTVSMALL)

        bool flag = false;
        iterations[faceI][0] = 0;

        for (int iterI=0; iterI<maxIter_; iterI++)
        {
            iterations[faceI][0] = iterI + 1;

            IntegrandFunc eddyViscosity =
                eddyViscosity_->value(sampler(), faceI, sqrt(tau), nuw[faceI]);

            scalar nuwI = nuw[faceI];

            IntegrandFunc integrand1 = [eddyViscosity, nuwI](const scalar y){
                return 1.0/(nuwI + eddyViscosity(y));
            };

            IntegrandFunc integrand2 = [eddyViscosity, nuwI](const scalar y){
                return y/(nuwI + eddyViscosity(y));
            };

            scalarList nutValues =
                eddyViscosity_->value(sampler(), faceI, y, sqrt(tau), nuw[faceI]);

            scalar integral1 = quad.integrate(integrand1, 0.0, sampler().h()[faceI], 1e-6);
            scalar integral2 = quad.integrate(integrand2, 0.0, sampler().h()[faceI], 1e-6);

            vector UFaceI(U[faceI][0], U[faceI][1], U[faceI][2]);

            scalar newTau =
                    sqr(magU[faceI]) + sqr(mag(sourceField[faceI])*integral2) -
                    2*(UFaceI & sourceField[faceI])*integral2;

            newTau  = sqrt(newTau)/(integral1 + VSMALL);

            scalar error = mag((tau - newTau) / tau);

            tau = newTau;
            

            if (iterI == maxIter_-1)
            {
                Pout << "havent converged" << nl;
                Pout << "Rank " << Pstream::myProcNo() << " Havent converged " << endl;
                WarningIn
                (
                    "Foam::ODEWallModelFvPatchScalarField::calcUTau()"
                )
                    << "tau_w did not converge to desired tolerance "
                    << eps_ << ". Error value: " << error << nl;
            }

            if (error < eps_)
            {
                flag = true;
                if (debug > 1)
                {
                    Info<< "tau_w converged after " << iterI + 1
                        << " iterations." << nl;
                }
                break;
            }

        }
        if (!flag)
        {
            Pout << "havent convergeed" << nl;
        }
            scalar tau_orig = tau;
            tau = old_tau;
            scalar tau_alt = 0;
            scalar up_alt = 0;

            for (int iterI=0; iterI<maxIter_; iterI++)
            {
                iterations[faceI][2] = iterI + 1;

                // New stuff
                IntegrandFunc nut_func = [kappa, Aplus](const Foam::scalar yPlus)
                {
                    return kappa * yPlus * Foam::sqr(1 - Foam::exp(-yPlus / Aplus));
                };

                IntegrandFunc integrand = [nut_func](const Foam::scalar yPlus)
                {
                    return 1.0 / (1.0 + nut_func(yPlus));
                };

                up_alt = quad.integrate(integrand, 0.0, sampler().h()[faceI] * sqrt(tau) / nuw[faceI], 1e-6);
                tau_alt = sqr(magU[faceI] / up_alt);

                scalar error = mag(tau - tau_alt)/tau;

                tau = tau_alt;
                
                if (error < eps_)
                {
                    break;
                }
            }

            tau = sqr(utau_exp);
            for (int iterI=0; iterI<maxIter_; iterI++)
            {
                iterations[faceI][1] = iterI + 1;

                IntegrandFunc eddyViscosity =
                    eddyViscosity_->value(sampler(), faceI, sqrt(tau), nuw[faceI]);

                scalar nuwI = nuw[faceI];

                IntegrandFunc integrand1 = [eddyViscosity, nuwI](const scalar y){
                    return 1.0/(nuwI + eddyViscosity(y));
                };

                IntegrandFunc integrand2 = [eddyViscosity, nuwI](const scalar y){
                    return y/(nuwI + eddyViscosity(y));
                };

                scalarList nutValues =
                    eddyViscosity_->value(sampler(), faceI, y, sqrt(tau), nuw[faceI]);

                scalar integral1 = quad.integrate(integrand1, 0.0, sampler().h()[faceI], 1e-6);
                scalar integral2 = quad.integrate(integrand2, 0.0, sampler().h()[faceI], 1e-6);

                vector UFaceI(U[faceI][0], U[faceI][1], U[faceI][2]);

                scalar newTau =
                        sqr(magU[faceI]) + sqr(mag(sourceField[faceI])*integral2) -
                        2*(UFaceI & sourceField[faceI])*integral2;

                newTau  = sqrt(newTau)/(integral1 + VSMALL);

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
                    WarningIn
                    (
                        "Foam::ODEWallModelFvPatchScalarField::calcUTau()"
                    )
                        << "tau_w did not converge to desired tolerance "
                        << eps_ << ". Error value: " << error << nl;
                }
            }

            //Pout << magU[faceI] << " " << sqrt(old_tau) << " " << sqrt(tau_alt) << " " << sqrt(tau_orig) << " " << utau_exp << " " << iterations[faceI][0] << " " << iterations[faceI][1] << " " << iterations[faceI][2] <<" " << (tau - sqr(utau_exp)) / tau * 100 << nl;

            uTau[faceI] = max(0.0, sqrt(tau_orig));

            wallModellingInfo[faceI][0] = tau_orig;
            wallModellingInfo[faceI][1] = sqr(utau_exp);
            wallModellingInfo[faceI][2] = tau_alt;

//        tau =  utau_exp;

    }

    if (patchSize)
    {
//        Pout << "Rrank " << Pstream::myProcNo() << " " << max(iterations) << endl;
    }

    //Pout << "iterations" << gMax(iterations) << gMin(iterations) << nl;



    // Grab global uTau field
    volScalarField & uTauField =
        const_cast<volScalarField &>
        (
            db().lookupObject<volScalarField>("uTauPredicted")
        );

    // Assign computed uTau to the boundary field of the global field
    uTauField.boundaryFieldRef()[patch().index()] == uTau;

    iterationsField.boundaryFieldRef()[patchi] == iterations;

//    Pout << gMax(iterationsField.boundaryFieldRef()[patchi]) << gMin(iterations) << nl;
    wallModellingInfoField.boundaryFieldRef()[patchi] == wallModellingInfo;

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
    sampler_(nullptr),
    meshes_(patch().size()),
    maxIter_(10),
    eps_(1e-3),
    nMeshY_(30)
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o1) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    createMeshes();
}

Foam::ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& orig,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallModelFvPatchScalarField(orig, p, iF, mapper),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    eddyViscosity_(orig.eddyViscosity_.clone()),
#else
    eddyViscosity_(orig.eddyViscosity_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler())),
    meshes_(orig.meshes_),
    maxIter_(orig.maxIter_),
    eps_(orig.eps_),
    nMeshY_(orig.nMeshY_)
{
    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o2) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    eddyViscosity_->addFieldsToSampler(sampler());
}


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
    sampler_
    (
        new SingleCellSampler
        (
            p,
            averagingTime(),
            dict.lookupOrDefault<word>("interpolationType", "cell"),
            dict.lookupOrDefault<word>("sampler", "Tree"),
            dict.lookupOrDefault<word>("lengthScaleType", "CubeRootVol"),
            dict.lookupOrDefault<bool>("hIsIndex", false)
        )
    ),
    meshes_(patch().size()),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 10)),
    eps_(dict.lookupOrDefault<scalar>("eps", 1e-3)),
    nMeshY_(dict.lookupOrDefault<label>("nMeshY", 30))

{
    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o3) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }

    createMeshes();
    eddyViscosity_->addFieldsToSampler(sampler());

    if (!db().found("solverIterations"))
    {
        db().store
        (
            new volVectorField
            (
                IOobject
                (
                    "solverIterations",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                vector(0,0,0),
                dimless
            )
        );
    }

    if (!db().found("wallModellingInfo"))
    {
        db().store
        (
            new volTensorField
            (
                IOobject
                (
                    "wallModellingInfo",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                tensor(0,0,0,0,0,0,0,0,0),
                dimless
            )
        );
    }
}


#ifdef FOAM_FVPATCHFIELD_NO_COPY
#else
Foam::ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& orig
)
:
    wallModelFvPatchScalarField(orig),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    eddyViscosity_(orig.eddyViscosity_.clone()),
#else
    eddyViscosity_(orig.eddyViscosity_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler())),
    meshes_(orig.meshes_),
    maxIter_(orig.maxIter_),
    eps_(orig.eps_),
    nMeshY_(orig.nMeshY_)

{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o4) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    eddyViscosity_->addFieldsToSampler(sampler());
}
#endif


Foam::ODEWallModelFvPatchScalarField::
ODEWallModelFvPatchScalarField
(
    const ODEWallModelFvPatchScalarField& orig,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallModelFvPatchScalarField(orig, iF),
#ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
    eddyViscosity_(orig.eddyViscosity_.clone()),
#else
    eddyViscosity_(orig.eddyViscosity_, false),
#endif
    sampler_(new SingleCellSampler(orig.sampler_())),
    meshes_(orig.meshes_),
    maxIter_(orig.maxIter_),
    eps_(orig.eps_),
    nMeshY_(orig.nMeshY_)
{

    if (debug)
    {
        Info<< "Constructing ODEWallModelfvPatchScalarField (o5) "
            << "from copy and DimensionedField for patch " << patch().name()
            << nl;
    }
    eddyViscosity_->addFieldsToSampler(sampler());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ODEWallModelFvPatchScalarField::write(Ostream& os) const
{
    wallModelFvPatchScalarField::write(os);
}


void Foam::ODEWallModelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    sampler().recomputeFields();
    sampler().sample();

    wallModelFvPatchScalarField::updateCoeffs();
}

// ************************************************************************* //
