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
    Foam::CrawlingCellFinder

@brief
    Class for searching for sampling cells by crawling through opposite faces.

    Starting from the wall-adjacent cell, the finder repeatedly identifies the
    face opposite to the incoming face and crosses it to the next cell. This is
    intended for meshes with layers of hexahedral or prismatic cells at the
    wall, where following opposite faces approximates the wall-normal
    direction. For other cells, no opposite face exists, and the crawl stops
    with a warning.

    The sampling height \f$h\f$ can be interpreted either as a physical
    distance from the wall face or as a consecutive off-wall cell index. In
    distance mode, crawling stops when the next opposing face is farther from
    the wall face than \f$h\f$, and the current cell is used. In index mode,
    \f$h\f$ is rounded to the requested layer index. For distances, the
    NormalRayCellFinder is the more general choice.

    When the opposing face lies on a processor boundary, the crawl is handed
    over to the neighbouring processor and continued there. The found sampling
    cells thus do not depend on the decomposition of the mesh. The crawl is
    stopped at all other boundaries, including processor boundaries
    originating from cyclic patches, and the last valid cell is used.
    Processor boundaries are also treated as other boundaries, when they are
    present outside a parallel run, e.g. in reconstructPar.

    For multi-cell sampling, every crawled cell is retained up to the requested
    distance or layer index. This search is restricted to the processor
    holding the wall face. If the crawl reaches a boundary patch before the
    requested target is reached, the last valid cells are retained. The
    optional \c excludeWallAdjacent argument removes the wall-adjacent cell
    only when at least one further cell remains.

    Invalid distances or indices fall back to the wall-adjacent cell.

    Used by default when \f$h\f$ is given as a cell index, see Sampler.

Contributors/Copyright:
    2019-2026 Timofey Mukha


SourceFiles
    CrawlingCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef CrawlingCellFinder_H
#define CrawlingCellFinder_H

#include "CellFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CrawlingCellFinder Declaration
\*---------------------------------------------------------------------------*/

class CrawlingCellFinder : public CellFinder
{

protected:

    // Protected Member Functions

        //- Crawl away from the wall within this processor, see CellFinder
        label advance
        (
            samplingLocation & location,
            const bool hIsIndex
        ) const override;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("CrawlingCellFinder");
#endif

    // Constructors

        //- Construct from patch
        CrawlingCellFinder
        (
            const fvPatch &
        );

        //- Construct from type and patch
        CrawlingCellFinder
        (
            const word & CrawlingCellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        CrawlingCellFinder
        (
            const CrawlingCellFinder &
        ) = default;

        //- Destructor
        virtual ~CrawlingCellFinder(){}

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
            const bool hIsIndex,
            const bool excludeWallAdjacent
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
