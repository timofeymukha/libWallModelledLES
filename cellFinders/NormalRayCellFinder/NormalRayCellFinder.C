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

#include "NormalRayCellFinder.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "codeRules.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(NormalRayCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, NormalRayCellFinder, Patch);
}
#endif

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::NormalRayCellFinder::intersectFace
(
    const label faceI,
    const samplingLocation & location,
    scalar & t,
    bool & leaving
) const
{
    const face & f = mesh_.faces()[faceI];
    const pointField & points = mesh_.points();
    const point & faceCentre = mesh_.faceCentres()[faceI];

    // The direction of the segment, into the domain
    const vector direction = -location.wallNormal;

    // Tolerance for the barycentric coordinates, to also accept segments
    // through the edges of the triangles
    const scalar tolerance = 1e-8;

    // Intersect with the triangles formed by the edges and the face centre
    forAll(f, i)
    {
        const vector e1 = points[f[i]] - faceCentre;
        const vector e2 = points[f.nextLabel(i)] - faceCentre;

        const vector p = direction ^ e2;
        const scalar det = e1 & p;

        // The segment is parallel to the triangle
        if (mag(det) <= SMALL*mag(e1)*mag(e2))
        {
            continue;
        }

        const vector s = location.wallCentre - faceCentre;
        const scalar u = (s & p)/det;

        if (u < -tolerance || u > 1 + tolerance)
        {
            continue;
        }

        const vector q = s ^ e1;
        const scalar v = (direction & q)/det;

        if (v < -tolerance || u + v > 1 + tolerance)
        {
            continue;
        }

        t = (e2 & q)/det;

        // The normal of the triangle, e1 ^ e2, points out of the owner, and
        // direction & (e1 ^ e2) = -det
        leaving = (det < 0) == (mesh_.faceOwner()[faceI] == location.cell);

        return true;
    }

    return false;
}


Foam::label Foam::NormalRayCellFinder::exitFacePlanes
(
    const samplingLocation & location,
    scalar & tIn,
    scalar & tOut
) const
{
    const cell & cellFaces = mesh_.cells()[location.cell];
    const vectorField & faceCentres = mesh_.faceCentres();
    const vectorField & faceAreas = mesh_.faceAreas();
    const labelList & owner = mesh_.faceOwner();

    // The direction of the segment, into the domain
    const vector direction = -location.wallNormal;

    tIn = -GREAT;
    tOut = GREAT;
    label exitFaceI = -1;

    // Clip the segment with the planes of the faces. The segment enters
    // through the faces with an outward normal against the direction, and
    // leaves through the ones with an outward normal along it.
    forAll(cellFaces, i)
    {
        const label faceI = cellFaces[i];

        const vector S =
            owner[faceI] == location.cell ? faceAreas[faceI] : -faceAreas[faceI];

        const scalar directionS = direction & S;

        // The segment is parallel to the face
        if (mag(directionS) <= SMALL*mag(S))
        {
            continue;
        }

        const scalar t =
            ((faceCentres[faceI] - location.wallCentre) & S)/directionS;

        if (directionS > 0)
        {
            if (t < tOut)
            {
                tOut = t;
                exitFaceI = faceI;
            }
        }
        else
        {
            tIn = max(tIn, t);
        }
    }

    return exitFaceI;
}


Foam::label Foam::NormalRayCellFinder::exitFace
(
    const samplingLocation & location,
    scalar & tIn,
    scalar & tOut
) const
{
    const cell & cellFaces = mesh_.cells()[location.cell];

    bool leaving;

    // The segment enters the cell through the face of the location
    if (!intersectFace(location.face, location, tIn, leaving))
    {
        return exitFacePlanes(location, tIn, tOut);
    }

    // Accept exits slightly before the entry, e.g. if the segment passes
    // through an edge of the cell
    const scalar tolerance = 1e-8*cbrt(mesh_.cellVolumes()[location.cell]);

    tOut = GREAT;
    label exitFaceI = -1;

    forAll(cellFaces, i)
    {
        const label faceI = cellFaces[i];
        scalar t;

        if
        (
            faceI != location.face
         && intersectFace(faceI, location, t, leaving)
         && leaving
         && t > tIn - tolerance
         && t < tOut
        )
        {
            tOut = t;
            exitFaceI = faceI;
        }
    }

    // Fall back to the planes of the faces, e.g. if the segment passes
    // between the triangles due to round-off
    if (exitFaceI == -1)
    {
        return exitFacePlanes(location, tIn, tOut);
    }

    return exitFaceI;
}


Foam::label Foam::NormalRayCellFinder::trace
(
    samplingLocation & location,
    DynamicList<label> & crossedCells,
    const bool parallel
) const
{
    const labelList & owner = mesh_.faceOwner();
    const labelList & neighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh();

    // Guard against cycling in degenerate cells
    const label maxLayers = 100000;

    while (location.layer < maxLayers)
    {
        crossedCells.append(location.cell);

        scalar tIn;
        scalar tOut;
        const label faceI = exitFace(location, tIn, tOut);

        if (faceI == -1)
        {
            Warning
                << "NormalRayCellFinder: Could not find the face through "
                << "which the wall-normal line of face "
                << location.originFace << " on patch " << patch().name()
                << " leaves cell " << location.cell << nl
                << "Will stop tracing and use this cell." << nl;

            return -1;
        }

        // Points closer to the exit face than the tolerance, relative to
        // the length of the segment within the cell, are on the face
        tIn = max(tIn, scalar(0));
        const scalar tolerance = 1e-8*(tOut - tIn);

        // The point at distance h is in this cell
        if (tOut > location.h + tolerance)
        {
            return -1;
        }

        if (mesh_.isInternalFace(faceI))
        {
            location.cell =
                owner[faceI] == location.cell
              ? neighbour[faceI]
              : owner[faceI];
            location.face = faceI;
            location.layer++;

            continue;
        }

        const polyPatch & exitPatch =
            boundaryMesh[boundaryMesh.whichPatch(faceI)];

        // Processor patches can also be present without a parallel run,
        // e.g. in reconstructPar, then they are treated as other patches
        if
        (
            parallel
         && Pstream::parRun()
         && isA<processorPolyPatch>(exitPatch)
         && !isA<processorCyclicPolyPatch>(exitPatch)
        )
        {
            location.cell = -1;
            location.face = faceI - exitPatch.start();
            location.layer++;

            return refCast<const processorPolyPatch>
            (
                exitPatch
            ).neighbProcNo();
        }

        // The segment leaves the domain, unless h is on the boundary,
        // sample in the middle of the part of the segment in the last cell
        if (tOut < location.h - tolerance)
        {
            const scalar h = 0.5*(tIn + tOut);

            Warning
                << "NormalRayCellFinder: The wall-normal line of face "
                << location.originFace << " on patch " << patch().name()
                << " leaves the domain through patch " << exitPatch.name()
                << " at distance " << tOut << ", before reaching the "
                << "requested distance " << location.h << nl
                << "Will stop tracing and use the last cell " << location.cell
                << " with the distance " << h << nl;

            location.h = h;
        }

        return -1;
    }

    Warning
        << "NormalRayCellFinder: Stopped tracing along the wall-normal line "
        << "of face " << location.originFace << " on patch " << patch().name()
        << " after " << maxLayers << " cells." << nl
        << "Will use the last cell " << location.cell << nl;

    return -1;
}


Foam::label Foam::NormalRayCellFinder::advance
(
    samplingLocation & location,
    const bool hIsIndex
) const
{
    DynamicList<label> crossedCells;
    return trace(location, crossedCells, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NormalRayCellFinder::NormalRayCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "NormalRayCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::NormalRayCellFinder::NormalRayCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::NormalRayCellFinder::findCells
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
            << "NormalRayCellFinder: hIsIndex is not supported by the "
            << "NormalRay sampler. Please use the Crawling sampler or "
            << "provide h as distances."
            <<  abort(FatalError);
    }

    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();

    // Locations, for which the sampling cell has been found
    DynamicList<samplingLocation> found(patch().size());

    // Locations, for which the trace continues on this processor
    DynamicList<samplingLocation> active(patch().size());

    forAll(patch(), patchFaceI)
    {
        const samplingLocation location =
            wallLocation(patchFaceI, h[patchFaceI], faceNormals);

        if (h[patchFaceI] < 0)
        {
            Warning
                << "NormalRayCellFinder: " << h[patchFaceI]
                << " is negative and thus not a valid distance. "
                << "Will fall back to wall-adjacent cell for face "
                <<  patchFaceI << " on patch " << patch().name() << nl;
            found.append(location);
        }
        else if (h[patchFaceI] == 0)
        {
            found.append(location);
        }
        else
        {
            active.append(location);
        }
    }

    search(active, found, hIsIndex);

    locations.transfer(found);
    createMap(locations, map);
}


void Foam::NormalRayCellFinder::findCellIndices
(
    labelListList & indexList,
    const scalarField & h,
    const bool excludeWallAdjacent
) const
{
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();

    DynamicList<label> crossedCells;

    forAll(patch(), patchFaceI)
    {
        samplingLocation location =
            wallLocation(patchFaceI, h[patchFaceI], faceNormals);

        if (h[patchFaceI] < 0)
        {
            Warning
                << "NormalRayCellFinder: " << h[patchFaceI]
                << " is negative and thus not a valid distance. "
                << "Will fall back to wall-adjacent cell for face "
                <<  patchFaceI << " on patch " << patch().name() << nl;
            indexList[patchFaceI] = labelList(1, location.cell);
            continue;
        }
        else if (h[patchFaceI] == 0)
        {
            indexList[patchFaceI] = labelList(1, location.cell);
            continue;
        }

        crossedCells.clear();
        trace(location, crossedCells, false);

        // Remove wall adjacent cell from list if applicable
        if (excludeWallAdjacent && crossedCells.size() > 1)
        {
            indexList[patchFaceI] = SubList<label>
            (
                crossedCells,
                crossedCells.size() - 1,
                1
            );
        }
        else
        {
            indexList[patchFaceI] = crossedCells;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
