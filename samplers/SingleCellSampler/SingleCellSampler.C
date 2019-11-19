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

#include "SingleCellSampler.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "SampledVelocityField.H"
#include "SampledWallGradUField.H"
#include "codeRules.H"
#include "patchDistMethod.H"
#include "scalarListIOList.H"
#include "CrawlingCellFinder.H"
#include "TreeCellFinder.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SingleCellSampler, 0);
    addToRunTimeSelectionTable(Sampler, SingleCellSampler, SamplerRTSTable);
}
#endif

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::SingleCellSampler::createIndexList()
{
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));

    
    scalarField hPatch = h.boundaryField()[patchIndex];

    if (hIsIndex_ && interpolationType_ != "cell")
    {
        Warning
            << "SingleCellSampler: hIsIndex is set to true, there is no sense "
            << "to interpolate within the cell." << nl
            << "Will fall back to the 'cell' interpolation type, i.e. use the "
            << "cell-centred values."
            << nl;

        interpolationType_ = "cell";
    }

    if (cellFinderType() == "Crawling")
    {
        CrawlingCellFinder cellFinder(patch());
        cellFinder.findCellIndices(indexList_, hPatch, hIsIndex_);
    }
    else if (cellFinderType() == word("Tree"))
    {
        if (hIsIndex_)
        {
            FatalErrorInFunction
                << "SingleCellSampler: hIsIndex is not supported by the Tree "
                << "sampler. Please use the Crawling sampler or provide h as "
                << "distances."
                <<  abort(FatalError);
        }
        TreeCellFinder cellFinder(patch());
        cellFinder.findCellIndices(indexList_, hPatch);
    }
    else
    {
        FatalErrorInFunction
            << "SingleCellSampler: cellFinderType should be either Tree or "
            << "Crawling. Current input is " << cellFinderType()
            <<  abort(FatalError);
    }


    const vectorField & patchFaceCentres = patch().Cf();
    const volVectorField & C = mesh_.C();


    // if the h field is an index, compute h_ as distance from the sampling
    // cell centers.
    // Same if h is distance, but we do not interpolate within the cell
    // Otherwise assign h_ to h.
    if (hIsIndex_ || (interpolationType() == "cell"))
    {
        forAll(patch(), faceI)
        {
            h_[faceI] = mag(C[indexList_[faceI]] - patchFaceCentres[faceI]);
        }
    }
    else
    {
        h_ = hPatch;
    }

    if (debug)
    {
        Info << "SingleCellSampler: Done" << nl;
    }
    
    // If the h field holds the distance, reassign the real distance used
    if (!hIsIndex())
    {
#ifdef FOAM_NEW_GEOMFIELD_RULES
        h.boundaryFieldRef()[patch().index()]
#else
        h.boundaryField()[patch().index()]
#endif
        ==
            h_;
    }

    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            mesh_.lookupObject<volScalarField> ("samplingCells")
        );
    
    forAll(indexList_, i)
    {
        samplingCells[indexList_[i]] = patchIndex; 
    }
}


void Foam::SingleCellSampler::createLengthList()
{
    // Cell volumes
    const scalarField & V = mesh_.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = pow(V[indexList_[i]], 1.0/3.0);
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SingleCellSampler::SingleCellSampler
(
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    bool hIsIndex
)
:
    Sampler(p, averagingTime, interpolationType, cellFinderType, hIsIndex),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0)
{
    createIndexList();
    createLengthList();
    
    addField
    (
            new SampledVelocityField(patch_, interpolationType_)
    );
    
    addField
    (
            new SampledWallGradUField(patch_, interpolationType_)
    );
}


Foam::SingleCellSampler::SingleCellSampler
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    bool hIsIndex
)
:
    SingleCellSampler
    (
        p,
        averagingTime,
        interpolationType,
        cellFinderType,
        hIsIndex
    )
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::SingleCellSampler::~SingleCellSampler()
{
}


void Foam::SingleCellSampler::sample() const
{
    // Ensure this processor has part of the patch
    if (!patch().size())
    {
        return;
    }    
    
    // Weight for time-averaging, default to 1 i.e no averaging.
    scalar eps = 1;
    if (averagingTime_ > mesh_.time().deltaTValue())
    {
        eps = mesh_.time().deltaTValue()/averagingTime_;
    }

    forAll(sampledFields_, fieldI)
    {

        scalarListList sampledList(patch().size());
        sampledFields_[fieldI].sample(sampledList, indexList(), h_);
        
        scalarListIOList & storedValues = const_cast<scalarListIOList & >
        (
            db().lookupObject<scalarListIOList>(sampledFields_[fieldI].name())
        );

        forAll(storedValues, i)
        {
            forAll(storedValues[i], j)
            {
                storedValues[i][j] = eps*sampledList[i][j] +
                                     (1 - eps)*storedValues[i][j];
            }
        }

    }
}


void Foam::SingleCellSampler::addField(SampledField * field)
{
    Sampler::addField(field);
    field->registerFields(indexList());
    field->setInterpolator(interpolationType_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
