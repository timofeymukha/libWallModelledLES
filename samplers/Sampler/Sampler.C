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
#include "codeRules.H"
#include "patchDistMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(Sampler, 0);
    defineRunTimeSelectionTable(Sampler, SamplerRTSTable);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::Sampler> Foam::Sampler::New 
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
{
    auto cstrIter =
    SamplerRTSTableConstructorTablePtr_->find(samplerName);

    if (cstrIter == SamplerRTSTableConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Sampler::New(const word&, const fvPatch & p, scalar averagingTime"
            
        )   << "Unknown Sampler type "
            << samplerName << nl << nl
            << "Valid Sampler types are :" << nl
            << SamplerRTSTableConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()
    (
        samplerName,
        p,
        averagingTime,
        interpolationType,
        cellFinderType,
        lengthScaleType,
        hIsIndex,
        excludeWallAdjacent
    );
}  

Foam::autoPtr<Foam::Sampler> Foam::Sampler::New 
(
    const dictionary & dict,
    const fvPatch & p
)
{
    word samplerName = dict.lookupOrDefault<word>("type", "SingleCellSampler");
    scalar averagingTime = dict.lookupOrDefault<scalar>("averagingTime", 0.0);
    word interpolationType =
        dict.lookupOrDefault<word>("interpolationType", "cell");
    word cellFinderType =
        dict.lookupOrDefault<word>("sampler", "Tree");
    word lengthScaleType =
        dict.lookupOrDefault<word>("lengthScale", "CubeRootVol");
    bool hIsIndex =
        dict.lookupOrDefault<bool>("hIsIndex", false);
    bool excludeWallAdjacent =
        dict.lookupOrDefault<bool>("excludeWallAdjacent", false);

    return Foam::Sampler::New
    (
        samplerName,
        p,
        averagingTime,
        interpolationType,
        cellFinderType,
        lengthScaleType,
        hIsIndex,
        excludeWallAdjacent
    );
}  


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Sampler::createFields()
{
    bool foundhSampler = mesh_.foundObject<volScalarField>("hSampler");
    word hName;

    // Grab h for the current patch
    if (foundhSampler)
    {
        hName = "hSampler";
    }
    else
    {
        hName = "h";
    }

    const volScalarField & h = mesh_.lookupObject<volScalarField> (hName);

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
                dimensionedScalar("samplingCells", dimless, -1),
                h.boundaryField().types()
            )
        );
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Sampler::Sampler
(
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
:
    patch_(p),
    averagingTime_(averagingTime),
    mesh_(patch_.boundaryMesh().mesh()),
    sampledFields_(0),
    interpolationType_(interpolationType),
    cellFinderType_(cellFinderType),
    lengthScaleType_(lengthScaleType),
    hIsIndex_(hIsIndex),
    excludeWallAdjacent_(excludeWallAdjacent)
{
    if (debug)
    {
        Info << "Sampler: Constructing from patch and avrg time" << nl;
    }

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

    if
    (
        !mesh_.subRegistry
        (
            "wallModelSampling"
        ).foundObject<objectRegistry>(patch_.name())
    )
    {
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
    }

    createFields();
}

Foam::Sampler::Sampler
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
:
    Sampler
    (
        p,
        averagingTime,
        interpolationType,
        cellFinderType,
        lengthScaleType,
        hIsIndex,
        excludeWallAdjacent
    )
{
}

Foam::Sampler::Sampler(const Sampler & copy)
:
    patch_(copy.patch_),
    averagingTime_(copy.averagingTime_),
    mesh_(copy.mesh_),
    sampledFields_(copy.sampledFields_),
    interpolationType_(copy.interpolationType_),
    cellFinderType_(copy.cellFinderType_),
    lengthScaleType_(copy.lengthScaleType_),
    hIsIndex_(copy.hIsIndex_),
    excludeWallAdjacent_(copy.excludeWallAdjacent_)
{
    if (debug)
    {
        Info << "Sampler: Running copy constructor" << nl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Sampler::~Sampler()
{
    if (debug)
    {
        Info << "Sampler: Running destructor" << nl;
    }
}

void Foam::Sampler::addField(SampledField * field)
{
    sampledFields_.setSize(sampledFields_.size() + 1);
    sampledFields_.set(sampledFields_.size() -1, field);
}


void Foam::Sampler::recomputeFields() const
{  
    forAll(sampledFields_, i)
    {
        sampledFields_[i].recompute();
    } 
}

void Foam::Sampler::write(Ostream & os) const
{
    os.writeKeyword("interpolationType")
        << interpolationType_ << token::END_STATEMENT << endl;
    os.writeKeyword("sampler")
        << cellFinderType_ << token::END_STATEMENT << endl;
    os.writeKeyword("hIsIndex")
        << hIsIndex_ << token::END_STATEMENT << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
