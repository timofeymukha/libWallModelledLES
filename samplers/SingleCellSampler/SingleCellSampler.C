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
#include "CellFinder.H"
#include "surfaceMesh.H"


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

    word hName;

    // Grab h for the current patch
    if (mesh_.foundObject<volScalarField>("hSampler"))
    {
        hName = "hSampler";
    }
    else
    {
        hName = "h";
    }

    // Grab h for the current patch
    volScalarField & h =
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> (hName));


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

    List<CellFinder::samplingLocation> locations;

    autoPtr<CellFinder> cellFinder
    (
        CellFinder::New(cellFinderType() + "CellFinder", patch())
    );
    cellFinder->findCells(locations, map_, hPatch, hIsIndex_);

    const vectorField & C = mesh_.cellCentres();
    const UList<label> & faceCells = patch().faceCells();
    const label myProc = Pstream::myProcNo();

    // Grab samplingCells field
    volScalarField & samplingCells =
        const_cast<volScalarField &>
        (
            mesh_.lookupObject<volScalarField> ("samplingCells")
        );

    sampledCells_.setSize(locations.size());
    sampledPoints_.setSize(locations.size());
    sampledWallCentres_.setSize(locations.size());

    // The distance from the wall used for each sampled cell
    scalarField sampledH(locations.size());

    indexList_ = -1;

    forAll(locations, i)
    {
        const CellFinder::samplingLocation & location = locations[i];
        const label cellI = location.cell;

        const bool wallAdjacent =
            location.originProc == myProc
         && cellI == faceCells[location.originFace];

        // if the h field is an index, compute h as distance from the sampling
        // cell centers.
        // Same if h is distance, but we do not interpolate within the cell
        // Otherwise use the prescribed h, excluding when we sample from the
        // wall-adjacent cell, since it can be used as fall-back when things
        // go wrong.
        if (hIsIndex_ || (interpolationType() == "cell") || wallAdjacent)
        {
            sampledH[i] = mag(C[cellI] - location.wallCentre);
        }
        else
        {
            sampledH[i] = location.h;
        }

        sampledCells_[i] = cellI;
        sampledWallCentres_[i] = location.wallCentre;
        sampledPoints_[i] =
            location.wallCentre - sampledH[i]*location.wallNormal;

        if (location.originProc == myProc)
        {
            indexList_[location.originFace] = cellI;
        }

        samplingCells[cellI] = patchIndex;
    }

    map_.distribute(sampledH);
    h_ = sampledH;

    if (debug)
    {
        Info << "SingleCellSampler: Done" << nl;
    }

    // If the global h field holds the distance, reassign the real distance used
    if (!hIsIndex())
    {
        h.boundaryFieldRef()[patch().index()] == h_;
    }
}


void Foam::SingleCellSampler::createLengthList(const word lengthScaleType)
{
    if (lengthScaleType == "CubeRootVol")
    {
        createLengthListCubeRootVol();
    }
    else if (lengthScaleType == "WallNormalDistance")
    {
        createLengthListWallNormalDistance();
    }
    else
    {
        FatalErrorInFunction
            << "SingleCellSampler: invalid length scale type, choose "
            << "CubeRootVol or WallNormalDistance. Current choice is "
            << lengthScaleType
            << abort(FatalError);
    }
}


void Foam::SingleCellSampler::createLengthListCubeRootVol()
{
    // Cell volumes
    const scalarField & V = mesh_.V();

    scalarField sampledLengths(sampledCells_.size());

    forAll(sampledLengths, i)
    {
        sampledLengths[i] = pow(V[sampledCells_[i]], 1.0/3.0);
    }

    map_.distribute(sampledLengths);
    lengthList_ = sampledLengths;
}


void Foam::SingleCellSampler::createLengthListWallNormalDistance()
{
    scalarField sampledLengths(sampledCells_.size());

    forAll(sampledLengths, i)
    {
        sampledLengths[i] =
            wallNormalDistance(sampledCells_[i], sampledWallCentres_[i]);
    }

    map_.distribute(sampledLengths);
    lengthList_ = sampledLengths;
}


Foam::scalar Foam::SingleCellSampler::wallNormalDistance
(
    const label index,
    const point & wallFaceCentre
) const
{
    const vectorField & faceCentres = mesh().faceCentres();
    const List<face> & faces = mesh().faces();
    const cell & cellI = mesh().cells()[index];

    // Find face with min distance from the wall face
    scalar minDist = GREAT;
    label minDistFace = 0;
    for (int j=0; j<cellI.nFaces(); j++)
    {
        const vector faceJ = faceCentres[cellI[j]];

        const scalar dist = mag(faceJ - wallFaceCentre);
        if (dist < minDist)
        {
            minDist = dist;
            minDistFace = cellI[j];
        }
    }

    label opposingFace =
        cellI.opposingFaceLabel(minDistFace, faces);
    vector opposingFaceCentre = faceCentres[opposingFace];
    vector minDistFaceCentre = faceCentres[minDistFace];

    return mag(opposingFaceCentre - minDistFaceCentre);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SingleCellSampler::SingleCellSampler
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
    Sampler(p, averagingTime, interpolationType, cellFinderType,
            lengthScaleType, hIsIndex, excludeWallAdjacent),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0),
    sampledCells_(),
    sampledPoints_(),
    sampledWallCentres_(),
    map_()
{
    if (!skipSamplingSetup())
    {
        createIndexList();
        createLengthList(lengthScaleType);
    }

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
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
:
    SingleCellSampler
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::SingleCellSampler::~SingleCellSampler()
{
}


void Foam::SingleCellSampler::sample() const
{
    if (skipSamplingSetup())
    {
        FatalErrorInFunction
            << "Sampling setup was skipped because "
            << "LIBWMLES_SKIP_SAMPLING_SETUP is enabled. "
            << "Unset it before running a solver that evaluates wall models."
            << exit(FatalError);
    }

    // No early return for empty patches, since this processor may hold
    // sampling cells for other processors.

    // Weight for time-averaging, default to 1 i.e no averaging.
    scalar eps = 1;
    if (averagingTime_ > mesh_.time().deltaTValue())
    {
        eps = mesh_.time().deltaTValue()/averagingTime_;
    }

    forAll(sampledFields_, fieldI)
    {
        const SampledField & field = sampledFields_[fieldI];

        scalarListList sampledList(patch().size());

        if (field.sampledFromCells())
        {
            // Interpolate at the sampled cells on this processor and send the
            // values to the processors holding the wall faces
            scalarListList values(sampledCells_.size());
            field.interpolate(values, sampledCells_, sampledPoints_);
            map_.distribute(values);
            field.postProcess(values);
            sampledList.transfer(values);
        }
        else
        {
            field.sample(sampledList, indexList(), h_);
        }

        scalarListIOList & storedValues = const_cast<scalarListIOList & >
        (
            db().lookupObject<scalarListIOList>(field.name())
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

    if (!skipSamplingSetup())
    {
        // The values are initialised from local cells, so use the
        // wall-adjacent cell for wall faces sampled from other processors.
        labelList cells(indexList());
        const UList<label> & faceCells = patch().faceCells();

        forAll(cells, i)
        {
            if (cells[i] == -1)
            {
                cells[i] = faceCells[i];
            }
        }

        field->registerFields(cells);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
