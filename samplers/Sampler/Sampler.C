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

#include "Sampler.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "SampledField.H"
#include "SampledVelocityField.H"
#include "SampledPGradField.H"
#include "SampledWallGradUField.H"
#include "codeRules.H"
#include "patchDistMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Sampler, 0);
    defineRunTimeSelectionTable(Sampler, PatchAndAveragingTime);
}
//
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::Sampler> Foam::Sampler::New 
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime
)
{
    PatchAndAveragingTimeConstructorTable::iterator cstrIter =
    PatchAndAveragingTimeConstructorTablePtr_->find(samplerName);

    if (cstrIter == WordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Sampler::New(const word&, const fvPatch & p, scalar averagingTime"
            
        )   << "Unknown Sampler type "
            << samplerName << nl << nl
            << "Valid Sampler types are :" << nl
            << WordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(samplerName, p, averagingTime);
}  

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Sampler::createFields()
{
      
    if (!mesh_.thisDb().found("h"))
    {
        if (debug)
        {
            Info<< "Sampler: Creating h field" << nl;
        }

        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));

    // Try to read field with distance from cell to patches
    
    if (debug)
    {
        Info<< "Sampler: Creating dist field" << nl;
    }

    if (!mesh_.thisDb().found("dist"))
    {
        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "dist",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("dist", dimLength,0),
                h.boundaryField().types()
            )
        );
    }
    volScalarField & dist = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField>("dist"));

    if (mag(max(dist.internalField()).value()) < VSMALL)
    {
        if (debug)
        {
            Info<< "Sampler: Computing dist field" << nl;
        }

        labelHashSet patchIDs(1);
        patchIDs.insert(patch().index());
        //Info<<patchIDs << nl;

        dictionary methodDict = dictionary();
        methodDict.lookupOrAddDefault(word("method"), word("meshWave"));

        if (debug)
        {
            Info<< "    Initializing patchDistanceMethod" << nl;
        }

        autoPtr<patchDistMethod>  pdm_
        (
            patchDistMethod::New
            (
                methodDict,
                mesh_,
                patchIDs
            )
        );

        if (debug)
        {
            Info<< "    Computing" << nl;
        }
        pdm_->correct(dist);
    }
    
   // Field that marks cells that are used for sampling
    if (!mesh_.thisDb().found("samplingCells"))
    {
        if (debug)
        {
            Info<< "Sampler: Creating samplingCells field" << nl;
        }

        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "samplingCells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("samplingCells", dimless,0),
                h.boundaryField().types()
            )
        );
    }
      
}


void Foam::Sampler::project(vectorField & field) const
{
    const tmp<vectorField> tfaceNormals = patch_.nf();
    const vectorField faceNormals = tfaceNormals();

    
    forAll(field, i)
    {    
        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(field[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        field[i] -= normal;        
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Sampler::Sampler
(
    const fvPatch & p,
    scalar averagingTime
)
:
    patch_(p),
    averagingTime_(averagingTime),
    mesh_(patch_.boundaryMesh().mesh()),
    sampledFields_(0)
{
    if (!mesh_.foundObject<objectRegistry>("wallModelSampling"))
    {
        objectRegistry * subObr = new objectRegistry
        (
            IOobject
            (
                "wallModelSampling",
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        subObr->store();
    }

    objectRegistry * subObr = new objectRegistry
    (
        IOobject
        (
            patch_.name(),
            mesh_.time().constant(),
            mesh_.subRegistry("wallModelSampling"),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    subObr->store();

    createFields();
}

Foam::Sampler::Sampler(const Sampler & copy)
:
    patch_(copy.patch_),
    averagingTime_(copy.averagingTime_),
    mesh_(copy.mesh_),
    sampledFields_(copy.sampledFields_.size())
{
    forAll(copy.sampledFields_, i)
    {
        sampledFields_[i] = copy.sampledFields_[i]->clone();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Sampler::~Sampler()
{
    forAll(sampledFields_, i)
    {
        delete sampledFields_[i];
    }
}

template<class Type>
void Foam::Sampler::listListToField
(
    const scalarListList & list,
    Field<Type> & field    
) const
{
    forAll(list, i)
    {
        Type element;
        forAll(list[i], j)
        {
            element[j] = list[i][j];
        }
        field[i] = element;
    }
}

void Foam::Sampler::addField(SampledField * field)
{
//    sampledFields_.append(field);
    sampledFields_.setSize(sampledFields_.size() + 1, field);
}


void Foam::Sampler::recomputeFields() const
{  
    forAll(sampledFields_, i)
    {
        sampledFields_[i]->recompute();
    } 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
