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

Class
    Foam::NormalRayCellFinder

@brief
    Class for finding sampling cells by tracing a ray along the wall normal.

    For each wall face, the sampling height \f$h\f$ is interpreted as a
    physical distance along the inward normal of the face. The finder follows
    the line segment from the face centre to the point at distance \f$h\f$,
    starting in the wall-adjacent cell. In each cell, the face through which
    the segment leaves the cell is found, and the trace continues in the cell
    on the other side of that face. The sampling cell is the cell containing
    the point at distance \f$h\f$. If this point lies exactly on a face
    between two cells, the cell further away from the wall is used.

    The faces are decomposed into triangles formed by their edges and their
    centres, the same as by OpenFOAM's face-centre decomposition of cells.
    The segment is intersected with these triangles, so warped faces and
    non-convex cells are handled consistently between neighbouring cells. No
    particular mesh structure is assumed, i.e. the finder works for arbitrary
    polyhedral cells, and not only for layered meshes. Interpreting \f$h\f$ as a
    consecutive off-wall cell index is not supported, use the
    CrawlingCellFinder for that.

    When the segment crosses a processor boundary, the trace is handed over
    to the neighbouring processor and continued there. The found sampling
    cells thus do not depend on the decomposition of the mesh. If the segment
    leaves the domain through any other boundary before reaching the
    distance \f$h\f$, including processor boundaries originating from cyclic
    patches, the last cell is used, and the sampling height is reduced to the
    middle of the part of the segment within that cell. Processor boundaries
    are also treated as other boundaries, when they are present outside a
    parallel run, e.g. in reconstructPar.

    For multi-cell sampling, all cells crossed by the segment are retained.
    This search is restricted to the processor holding the wall face. The
    optional \c excludeWallAdjacent argument removes the wall-adjacent cell
    only when at least one further cell remains.

    Distances \f$h \le 0\f$ fall back to the wall-adjacent cell.

    Used by default when \f$h\f$ is given as a distance, see Sampler.

Contributors/Copyright:
    2026 Timofey Mukha


SourceFiles
    NormalRayCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef NormalRayCellFinder_H
#define NormalRayCellFinder_H

#include "CellFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NormalRayCellFinder Declaration
\*---------------------------------------------------------------------------*/

class NormalRayCellFinder : public CellFinder
{

protected:

    // Protected Member Functions

        //- Intersect the wall-normal line of a location with a face of its
        //- cell. The face is decomposed into triangles formed by its edges
        //- and its centre. Returns whether the line intersects the face, the
        //- distance from the wall to the intersection, and whether the line
        //- leaves the cell there.
        bool intersectFace
        (
            const label faceI,
            const samplingLocation & location,
            scalar & t,
            bool & leaving
        ) const;

        //- As exitFace, but treating the faces as planes through their
        //- centres, which is exact for convex cells with flat faces
        label exitFacePlanes
        (
            const samplingLocation & location,
            scalar & tIn,
            scalar & tOut
        ) const;

        //- The part of the wall-normal segment of a location, which lies
        //- within its cell, given by the distances from the wall, at which
        //- the segment enters and leaves the cell. Returns the face through
        //- which the segment leaves the cell, or -1 if none is found.
        label exitFace
        (
            const samplingLocation & location,
            scalar & tIn,
            scalar & tOut
        ) const;

        //- Trace the ray along the wall normal within this processor,
        //- appending the crossed cells. Returns -1 if the trace is complete,
        //- otherwise the neighbouring processor to which it should be handed
        //- over. If not parallel, processor boundaries are treated as other
        //- boundaries.
        label trace
        (
            samplingLocation & location,
            DynamicList<label> & crossedCells,
            const bool parallel
        ) const;

        //- Trace the ray along the wall normal within this processor, see
        //- CellFinder
        label advance
        (
            samplingLocation & location,
            const bool hIsIndex
        ) const override;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("NormalRayCellFinder");
#endif

    // Constructors

        //- Construct from patch
        NormalRayCellFinder
        (
            const fvPatch &
        );

        //- Construct from type and patch
        NormalRayCellFinder
        (
            const word & NormalRayCellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        NormalRayCellFinder
        (
            const NormalRayCellFinder &
        ) = default;

        //- Destructor
        virtual ~NormalRayCellFinder(){}

    // Member functions

        //- Find a sampling cell for each face of the patch, possibly on
        //- another processor, see CellFinder
        void findCells
        (
            List<samplingLocation> & locations,
            mapDistributeBase & map,
            const scalarField & h,
            const bool hIsIndex
        ) const override;

        //- Find the sampling cells for multi-cell sampling for each face of
        //- the patch, searching on this processor only
        void findCellIndices
        (
            labelListList & indexList,
            const scalarField & h,
            const bool excludeWallAdjacent
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
