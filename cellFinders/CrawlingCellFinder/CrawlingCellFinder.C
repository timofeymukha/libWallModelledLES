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

#include "CrawlingCellFinder.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "codeRules.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CrawlingCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, CrawlingCellFinder, Patch);
}
#endif

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::CrawlingCellFinder::advance
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

    // For distances, an arbitrary big number
    const label nLayers = hIsIndex ? label(std::round(location.h)) : 100;

    while (location.layer < nLayers)
    {
        const label opposingFace =
            cells[location.cell].opposingFaceLabel(location.face, faces);

        if (opposingFace == -1)
        {
            Warning
                << "CrawlingCellFinder: Could not find opposing face "
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
                    << "CrawlingCellFinder: The opposing face for cell "
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CrawlingCellFinder::CrawlingCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "CrawlingCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::CrawlingCellFinder::CrawlingCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CrawlingCellFinder::findCells
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
                    << "CrawlingCellFinder: Could not round "
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
                    << "CrawlingCellFinder: " << h[patchFaceI]
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

    search(active, found, hIsIndex);

    locations.transfer(found);
    createMap(locations, map);
}


void Foam::CrawlingCellFinder::findCellIndices
(
    labelListList & indexList,
    const scalarField & h,
    const bool hIsIndex,
    const bool excludeWallAdjacent
) const
{
    const labelList & owner = mesh_.faceOwner();
    const labelList & neighbour = mesh_.faceNeighbour();
    const vectorField & faceCentres = mesh_.faceCentres();
    const List<cell> & cells = mesh().cells();

    const vectorField & patchFaceCentres = patch().Cf();
    const vectorField & C = mesh_.cellCentres();
    const UList<label> & faceCells = patch().faceCells();
    const List<face> & faces = mesh().faces();
    const label patchStart = patch().start();
    const polyBoundaryMesh & boundaryMesh = mesh().boundaryMesh();

    forAll(patch(), patchFaceI)
    {

        label nLayers = 0;

        if (hIsIndex)
        {
            nLayers = std::round(h[patchFaceI]);
            if (nLayers < 1)
            {
                Warning
                    << "CrawlingCellFinder: Could not round "
                    << h[patchFaceI] << " to a valid cell layer index. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
            else if (nLayers == 1)
            {
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
        }
        else
        {
            if (h[patchFaceI] < 0)
            {
                Warning
                    << "CrawlingCellFinder: " << h[patchFaceI]
                    << " is negative and thus not a valid distance. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
            else if (h[patchFaceI] == 0)
            {
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }

            nLayers = 100; // arbitrary big number
        }
        indexList[patchFaceI].setSize(nLayers, -1);

        label startCellIndex = faceCells[patchFaceI];
        cell startCell = cells[startCellIndex];
        label startFaceLabel = patchStart + patchFaceI;
        label nextCellIndex = -1;

        label layerCounter = 0;

        for (label layer=0; layer < nLayers; layer++)
        {
            layerCounter++;

            // Grab the next cell
            indexList[patchFaceI][layer] = startCellIndex;

            label opposingFace = 
                startCell.opposingFaceLabel(startFaceLabel, faces);

            if (opposingFace == -1)
            {
                Warning
                    << "CrawlingCellFinder: Could not find opposing face"
                    << "for cell " << startCellIndex << " with cell center "
                    << C[startCellIndex] << " and face " << startFaceLabel << nl
                    << "Will stop crawling and use the last valid cell " 
                    << "corresponding to index " << layer + 1
                    << nl;
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }
            // if we hit a patch
            else if (opposingFace >= mesh().nInternalFaces())
            {
                label opposingPatchInd = boundaryMesh.whichPatch(opposingFace);
                const fvPatch & opposingPatch =
                    mesh().boundary()[opposingPatchInd];
                word opposingPatchName = opposingPatch.name();

                // Recompute as the actual distance that will be used
                // due to abortion in crawling
                scalar distance =
                    mag(C[startCellIndex] - patchFaceCentres[patchFaceI]);

                Warning
                    << "CrawlingCellFinder: The opposing face for cell "
                    << startCellIndex << " with cell center "
                    << C[startCellIndex] << " and face " << startFaceLabel
                    << " belongs to patch " << opposingPatchName << nl 
                    << "Will stop crawling and use the last valid cell " 
                    << "corresponding to index " << layer + 1
                    << " and distance " << distance
                    << nl;
                    
                // Need this when hIsIndex
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }

            scalar distance =
                mag(faceCentres[opposingFace] - patchFaceCentres[patchFaceI]);

            // if the opposing face is above h, stop 
            if ((!hIsIndex) && (distance > h[patchFaceI]))
            {
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }

            // find index of the next cell
            if (owner[opposingFace] == startCellIndex)
            {
                nextCellIndex = neighbour[opposingFace];
            }
            else
            {
                nextCellIndex = owner[opposingFace];
            }

            startFaceLabel = opposingFace;
            startCellIndex = nextCellIndex;
            startCell = cells[startCellIndex];
        }

        if ((indexList[patchFaceI].size() != 1) && excludeWallAdjacent)
        {
            for(int j=0; j<indexList[patchFaceI].size()-1; j++)
            {
                indexList[patchFaceI][j] = indexList[patchFaceI][j+1];
            }
            indexList[patchFaceI].setSize(indexList[patchFaceI].size()-1);

        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
