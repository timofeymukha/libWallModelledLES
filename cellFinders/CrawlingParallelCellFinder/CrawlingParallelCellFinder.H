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
    Foam::CrawlingParallelCellFinder

@brief
    Class for finding sampling cells by crawling away from the wall, also
    across processor boundaries.

    The search is the same as for the CrawlingCellFinder, see its
    documentation. However, when the opposing face lies on a processor
    boundary, the crawl is handed over to the neighbouring processor and
    continued there. The found sampling cells thus do not depend on the
    decomposition of the mesh. The crawl is stopped at all other boundaries,
    including processor boundaries originating from cyclic patches, the same
    as in serial. Processor boundaries are also treated as other boundaries,
    when they are present outside a parallel run, e.g. in reconstructPar.

    Selected with the sampler keyword set to CrawlingParallel.

Contributors/Copyright:
    2026 Timofey Mukha


SourceFiles
    CrawlingParallelCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef CrawlingParallelCellFinder_H
#define CrawlingParallelCellFinder_H

#include "CellFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CrawlingParallelCellFinder Declaration
\*---------------------------------------------------------------------------*/

class CrawlingParallelCellFinder : public CellFinder
{

protected:

    // Protected Member Functions

        //- Crawl away from the wall within this processor. Returns -1 if
        //- the sampling cell has been found, otherwise the neighbouring
        //- processor to which the crawl should be handed over.
        label crawl(samplingLocation & location, const bool hIsIndex) const;

        //- For each processor, the index of the processor patch shared
        //- with it, or -1. Processor patches stemming from cyclics are
        //- not considered.
        labelList processorPatches() const;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("CrawlingParallelCellFinder");
#endif

    // Constructors

        //- Construct from patch
        CrawlingParallelCellFinder
        (
            const fvPatch &
        );

        //- Construct from type and patch
        CrawlingParallelCellFinder
        (
            const word & CrawlingParallelCellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        CrawlingParallelCellFinder
        (
            const CrawlingParallelCellFinder &
        ) = default;

        //- Destructor
        virtual ~CrawlingParallelCellFinder(){}

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
