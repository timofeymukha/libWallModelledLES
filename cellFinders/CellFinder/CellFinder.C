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

#include "CellFinder.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "codeRules.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CellFinder, 0);
    defineRunTimeSelectionTable(CellFinder, Patch);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::CellFinder> Foam::CellFinder::New 
(
    const word & CellFinderName,
    const fvPatch & p
)
{
    if (CellFinderName == "TreeCellFinder")
    {
        FatalErrorInFunction
            << "The Tree sampler has been removed. Remove the sampler "
            << "keyword, then the NormalRay sampler is used for distances, "
            << "which, like the Tree sampler, finds the cell at distance h "
            << "along the wall normal."
            << exit(FatalError);
    }

    auto cstrIter =
    PatchConstructorTablePtr_->find(CellFinderName);

    if (cstrIter == PatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "CellFinder::New(const word&, const fvPatch & p"
            
        )   << "Unknown CellFinder type "
            << CellFinderName << nl << nl
            << "Valid CellFinder types are :" << nl
            << PatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(CellFinderName, p);
}  

Foam::autoPtr<Foam::CellFinder> Foam::CellFinder::New 
(
    const dictionary & dict,
    const fvPatch & p
)
{
    word CellFinderName =
        dict.lookupOrDefault<word>("type", "SingleCellCellFinder");

    return Foam::CellFinder::New(CellFinderName, p);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::CellFinder::createFields() const
{
}


Foam::CellFinder::samplingLocation Foam::CellFinder::wallLocation
(
    const label patchFaceI,
    const scalar h,
    const vectorField & faceNormals
) const
{
    samplingLocation location;
    location.originProc = Pstream::myProcNo();
    location.originFace = patchFaceI;
    location.cell = patch().faceCells()[patchFaceI];
    location.face = patch().start() + patchFaceI;
    location.layer = 1;
    location.wallCentre = patch().Cf()[patchFaceI];
    location.wallNormal = faceNormals[patchFaceI];
    location.h = h;

    return location;
}


Foam::labelList Foam::CellFinder::processorPatches() const
{
    const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh();

    labelList patches(Pstream::nProcs(), -1);

    if (!Pstream::parRun())
    {
        return patches;
    }

    forAll(boundaryMesh, patchI)
    {
        const polyPatch & pp = boundaryMesh[patchI];

        if (isA<processorPolyPatch>(pp) && !isA<processorCyclicPolyPatch>(pp))
        {
            patches[refCast<const processorPolyPatch>(pp).neighbProcNo()] =
                patchI;
        }
    }

    return patches;
}


void Foam::CellFinder::search
(
    DynamicList<samplingLocation> & active,
    DynamicList<samplingLocation> & found,
    const bool hIsIndex
) const
{
    const labelList procPatches(processorPatches());
    const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh();
    const labelList & owner = mesh_.faceOwner();

    while (true)
    {
        // Searches to be handed over to the neighbouring processors
        List<DynamicList<samplingLocation>> handovers(Pstream::nProcs());

        forAll(active, i)
        {
            const label nbrProc = advance(active[i], hIsIndex);

            if (nbrProc == -1)
            {
                found.append(active[i]);
            }
            else
            {
                handovers[nbrProc].append(active[i]);
            }
        }

        active.clear();

        if (!Pstream::parRun())
        {
            break;
        }

        List<List<samplingLocation>> received;
        exchange(handovers, received);

        forAll(received, proci)
        {
            if (received[proci].size() && procPatches[proci] == -1)
            {
                FatalErrorInFunction
                    << "CellFinder: Received searches from processor "
                    << proci << ", which is not a neighbour."
                    << exit(FatalError);
            }

            // The face is entered from the neighbouring processor, continue
            // in the cell on this side
            forAll(received[proci], i)
            {
                samplingLocation & location = received[proci][i];

                location.face += boundaryMesh[procPatches[proci]].start();
                location.cell = owner[location.face];
                active.append(location);
            }
        }

        if (returnReduce(active.size(), sumOp<label>()) == 0)
        {
            break;
        }
    }
}


void Foam::CellFinder::createMap
(
    const UList<samplingLocation> & locations,
    mapDistributeBase & map
) const
{
    const label nProcs = Pstream::nProcs();

    // Per processor holding the wall faces, the locations to send there and
    // the corresponding wall faces
    List<DynamicList<label>> sendLocations(nProcs);
    List<DynamicList<label>> originFaces(nProcs);

    forAll(locations, i)
    {
        sendLocations[locations[i].originProc].append(i);
        originFaces[locations[i].originProc].append(locations[i].originFace);
    }

    labelListList constructMap;
    exchange(originFaces, constructMap);

    // Every wall face should receive exactly one sample
    labelList nSamples(patch().size(), 0);

    forAll(constructMap, proci)
    {
        forAll(constructMap[proci], i)
        {
            nSamples[constructMap[proci][i]]++;
        }
    }

    forAll(nSamples, faceI)
    {
        if (nSamples[faceI] != 1)
        {
            FatalErrorInFunction
                << "CellFinder: Face " << faceI << " on patch "
                << patch().name() << " has " << nSamples[faceI]
                << " sampling cells instead of exactly one."
                << exit(FatalError);
        }
    }

    labelListList subMap(nProcs);

    forAll(subMap, proci)
    {
        subMap[proci].transfer(sendLocations[proci]);
    }

    map = mapDistributeBase
    (
        patch().size(),
        std::move(subMap),
        std::move(constructMap)
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CellFinder::CellFinder
(
    const fvPatch & p
)
:
    patch_(p),
    mesh_(patch_.boundaryMesh().mesh())
{
    if (debug)
    {
        Info << "CellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::CellFinder::CellFinder
(
    const word & CellFinderName,
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "CellFinder: Constructing from name and patch" << nl;
    }
}

Foam::CellFinder::CellFinder(const CellFinder & copy)
:
    patch_(copy.patch_),
    mesh_(copy.mesh_)
{
    if (debug)
    {
        Info << "CellFinder: Running copy constructor" << nl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::CellFinder::~CellFinder()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
