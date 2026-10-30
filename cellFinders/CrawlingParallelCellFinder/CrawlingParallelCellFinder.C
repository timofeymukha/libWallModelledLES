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

#include "CrawlingParallelCellFinder.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "codeRules.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CrawlingParallelCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, CrawlingParallelCellFinder, Patch);
}
#endif

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::CrawlingParallelCellFinder::crawl
(
    samplingLocation & location,
    const bool hIsIndex
) const
{
    const labelList & owner = mesh_.faceOwner();
    const labelList & neighbour = mesh_.faceNeighbour();
    const vectorField & faceCentres = mesh_.faceCentres();
    const vectorField & C = mesh_.cellCentres();
    const List<cell> & cells = mesh_.cells();
    const List<face> & faces = mesh_.faces();
    const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh();

    // For distances, an arbitrary big number, as in the CrawlingCellFinder
    const label nLayers = hIsIndex ? label(std::round(location.h)) : 100;

    while (location.layer < nLayers)
    {
        const label opposingFace =
            cells[location.cell].opposingFaceLabel(location.face, faces);

        if (opposingFace == -1)
        {
            Warning
                << "CrawlingParallelCellFinder: Could not find opposing face "
                << "for cell " << location.cell << " with cell center "
                << C[location.cell] << " and face " << location.face << nl
                << "Will stop crawling and use the last valid cell "
                << "corresponding to index " << location.layer
                << nl;

            return -1;
        }

        const scalar distance =
            mag(faceCentres[opposingFace] - location.wallCentre);

        if (opposingFace >= mesh_.nInternalFaces())
        {
            const polyPatch & opposingPatch =
                boundaryMesh[boundaryMesh.whichPatch(opposingFace)];

            // Processor patches can also be present without a parallel run,
            // e.g. in reconstructPar, then they are treated as other patches
            if
            (
                Pstream::parRun()
             && isA<processorPolyPatch>(opposingPatch)
             && !isA<processorCyclicPolyPatch>(opposingPatch)
            )
            {
                // The opposing face is internal in the undecomposed mesh,
                // so check the distance as for any other internal face.
                if ((!hIsIndex) && (distance > location.h))
                {
                    return -1;
                }

                location.cell = -1;
                location.face = opposingFace - opposingPatch.start();
                location.layer++;

                return refCast<const processorPolyPatch>
                (
                    opposingPatch
                ).neighbProcNo();
            }

            // distance that will be used due to abortion in crawling
            const scalar cellDistance =
                mag(C[location.cell] - location.wallCentre);

            // If distance-based we might actually get the right h
            // in the last cell, so need to check for that
            // before issuing the warning
            if
            (
                hIsIndex
             || mag(cellDistance - location.h)/cellDistance > SMALL
            )
            {
                Warning
                    << "CrawlingParallelCellFinder: The opposing face for cell "
                    << location.cell << " with cell center "
                    << C[location.cell] << " and face " << location.face
                    << " belongs to patch " << opposingPatch.name() << nl
                    << "Will stop crawling and use the last valid cell "
                    << "corresponding to index " << location.layer
                    << " and cell centre at distance " << cellDistance
                    << ". Requested index/distance is " << location.h
                    << " which may or may not correspond to this cell."
                    << nl;
            }

            return -1;
        }

        // if the opposing face is above h, stop and grab the cell
        if ((!hIsIndex) && (distance > location.h))
        {
            return -1;
        }

        if (owner[opposingFace] == location.cell)
        {
            location.cell = neighbour[opposingFace];
        }
        else
        {
            location.cell = owner[opposingFace];
        }

        location.face = opposingFace;
        location.layer++;
    }

    return -1;
}


Foam::labelList Foam::CrawlingParallelCellFinder::processorPatches() const
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CrawlingParallelCellFinder::CrawlingParallelCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "CrawlingParallelCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::CrawlingParallelCellFinder::CrawlingParallelCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CrawlingParallelCellFinder::findCells
(
    List<samplingLocation> & locations,
    mapDistributeBase & map,
    const scalarField & h,
    const bool hIsIndex
) const
{
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();

    // Locations, for which the sampling cell has been found
    DynamicList<samplingLocation> found(patch().size());

    // Locations, for which the crawling continues on this processor
    DynamicList<samplingLocation> active(patch().size());

    forAll(patch(), patchFaceI)
    {
        const samplingLocation location =
            wallLocation(patchFaceI, h[patchFaceI], faceNormals);

        if (hIsIndex)
        {
            const label nLayers = std::round(h[patchFaceI]);

            if (nLayers < 1)
            {
                Warning
                    << "CrawlingParallelCellFinder: Could not round "
                    << h[patchFaceI] << " to a valid cell layer index. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                found.append(location);
                continue;
            }
            else if (nLayers == 1)
            {
                found.append(location);
                continue;
            }
        }
        else
        {
            if (h[patchFaceI] < 0)
            {
                Warning
                    << "CrawlingParallelCellFinder: " << h[patchFaceI]
                    << " is negative and thus not a valid distance. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                found.append(location);
                continue;
            }
            else if (h[patchFaceI] == 0)
            {
                found.append(location);
                continue;
            }
        }

        active.append(location);
    }

    const labelList procPatches(processorPatches());
    const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh();
    const labelList & owner = mesh_.faceOwner();

    while (true)
    {
        // Crawls to be handed over to the neighbouring processors
        List<DynamicList<samplingLocation>> handovers(Pstream::nProcs());

        forAll(active, i)
        {
            const label nbrProc = crawl(active[i], hIsIndex);

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
                    << "CrawlingParallelCellFinder: Received crawls from "
                    << "processor " << proci << ", which is not a neighbour."
                    << exit(FatalError);
            }

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

    locations.transfer(found);
    createMap(locations, map);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
