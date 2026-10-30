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

#include "TreeParallelCellFinder.H"
#include "TreeCellFinder.H"
#include "volFields.H"
#include "treeDataCell.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "codeRules.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(TreeParallelCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, TreeParallelCellFinder, Patch);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TreeParallelCellFinder::TreeParallelCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "TreeParallelCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::TreeParallelCellFinder::TreeParallelCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TreeParallelCellFinder::findCells
(
    List<samplingLocation> & locations,
    mapDistributeBase & map,
    const scalarField & h,
    const bool hIsIndex
) const
{
    if (hIsIndex)
    {
        FatalErrorInFunction
            << "TreeParallelCellFinder: hIsIndex is not supported by the "
            << "TreeParallel sampler. Please use the CrawlingParallel sampler "
            << "or provide h as distances."
            <<  abort(FatalError);
    }

    const label nProcs = Pstream::nProcs();
    const label myProc = Pstream::myProcNo();

    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();

    const scalar maxH = gMax(h);

    // Construct the octree over the candidate cells on this processor
    labelField candidates;
    autoPtr<indexedOctree<treeDataCell> > treePtr;

    if (maxH > 0)
    {
        TreeCellFinder treeFinder(patch());
        tmp<volScalarField> distField(treeFinder.distanceField());

        candidates = treeFinder.findCandidateCellLabels
        (
            distField().primitiveField(),
            scalarField(1, maxH)
        );

        // The tet base points are used when searching for the cell
        // containing a point. They are computed on demand with parallel
        // communication, so trigger it here on all processors, since only
        // some will receive points to search for.
        (void)mesh_.tetBasePtIs();

        if (candidates.size())
        {
            treeBoundBox boundBox(mesh_.bounds());
            Random rndGen(261782);
#ifdef FOAM_TREEBOUNDBOX_DOES_NOT_ACCEPT_RNG
            boundBox.extend(1e-4);
#else
            boundBox.extend(rndGen, 1e-4);
#endif
            boundBox.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            boundBox.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

            treePtr.reset
            (
                new indexedOctree<treeDataCell>
                (
                    treeDataCell
                    (
                        false,
                        mesh_,
                        candidates,
                        polyMesh::CELL_TETS
                    ),
                    boundBox,
                    8,
                    10,
                    3.0
                )
            );
        }
    }

    // The bounding boxes of the searchable cells on all processors
    List<boundBox> procBoundBoxes(nProcs, boundBox::invertedBox);

    if (treePtr)
    {
        procBoundBoxes[myProc] = treePtr->bb();
    }

    Pstream::allGatherList(procBoundBoxes);

    // Locations, for which the sampling cell has been found
    DynamicList<samplingLocation> found(patch().size());

    // Points to search for, sent to each processor
    List<DynamicList<samplingLocation>> queries(nProcs);
    boolList queried(patch().size(), false);

    forAll(patch(), patchFaceI)
    {
        const samplingLocation location =
            wallLocation(patchFaceI, h[patchFaceI], faceNormals);

        if (h[patchFaceI] < 0)
        {
            Warning
                << "TreeParallelCellFinder: " << h[patchFaceI]
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

        const point p =
            location.wallCentre - location.wallNormal*location.h;

        forAll(procBoundBoxes, proci)
        {
            if (procBoundBoxes[proci].contains(p))
            {
                queries[proci].append(location);
                queried[patchFaceI] = true;
            }
        }

        if (!queried[patchFaceI])
        {
            Warning
                << "TreeParallelCellFinder: the point " << h[patchFaceI]
                << " away from the wall is outside the domain. "
                << "Will fall back to wall-adjacent cell for face "
                <<  patchFaceI << " on patch " << patch().name() << nl;
            found.append(location);
        }
    }

    List<List<samplingLocation>> received;
    exchange(queries, received);

    // Search for the received points, -1 if not found
    List<DynamicList<label>> answers(nProcs);

    forAll(received, proci)
    {
        forAll(received[proci], i)
        {
            const samplingLocation & location = received[proci][i];

            label cellI = -1;

            if (treePtr)
            {
                const label index = treePtr->findInside
                (
                    location.wallCentre - location.wallNormal*location.h
                );

                if (index != -1)
                {
                    cellI = candidates[index];
                }
            }

            answers[proci].append(cellI);
        }
    }

    labelListList replies;
    exchange(answers, replies);

    // Select the lowest processor that found the point
    labelList selectedProc(patch().size(), -1);
    labelList selectedCell(patch().size(), -1);

    forAll(replies, proci)
    {
        forAll(replies[proci], i)
        {
            const label patchFaceI = queries[proci][i].originFace;

            if (replies[proci][i] != -1 && selectedProc[patchFaceI] == -1)
            {
                selectedProc[patchFaceI] = proci;
                selectedCell[patchFaceI] = replies[proci][i];
            }
        }
    }

    // Send the selection to the processors holding the sampling cells
    List<DynamicList<samplingLocation>> selections(nProcs);

    forAll(patch(), patchFaceI)
    {
        if (!queried[patchFaceI])
        {
            continue;
        }

        samplingLocation location =
            wallLocation(patchFaceI, h[patchFaceI], faceNormals);

        if (selectedProc[patchFaceI] == -1)
        {
            Warning
                << "TreeParallelCellFinder: no cell containing the point "
                << h[patchFaceI] << " away from the wall was found. "
                << "Will fall back to wall-adjacent cell for face "
                <<  patchFaceI << " on patch " << patch().name() << nl;
            found.append(location);
        }
        else
        {
            location.cell = selectedCell[patchFaceI];
            selections[selectedProc[patchFaceI]].append(location);
        }
    }

    List<List<samplingLocation>> selected;
    exchange(selections, selected);

    forAll(selected, proci)
    {
        found.append(selected[proci]);
    }

    locations.transfer(found);
    createMap(locations, map);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
