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
    Foam::TreeParallelCellFinder

@brief
    Class for finding the cells containing the points at the sampling
    distance from the wall, also across processor boundaries.

    For each wall face, the point \f$h\f$ away from the face centre along the
    inward normal is sent to all processors whose bounding box contains it,
    and the cell containing the point is searched for using an octree over
    the candidate cells, see the TreeCellFinder. If several processors find
    the point, e.g. because it lies on a processor boundary, the one with
    the lowest rank is selected. If no cell is found, the point is outside
    the domain and the wall-adjacent cell is used. The found sampling cells
    thus do not depend on the decomposition of the mesh.

    Note that the TreeCellFinder selects the cell with the nearest cell
    centre instead, which in rare cases can differ from the cell containing
    the point. Interpreting \f$h\f$ as a cell index is not supported.

    Selected with the sampler keyword set to TreeParallel.

Contributors/Copyright:
    2026 Timofey Mukha


SourceFiles
    TreeParallelCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef TreeParallelCellFinder_H
#define TreeParallelCellFinder_H

#include "CellFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TreeParallelCellFinder Declaration
\*---------------------------------------------------------------------------*/

class TreeParallelCellFinder : public CellFinder
{

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("TreeParallelCellFinder");
#endif

    // Constructors

        //- Construct from patch
        TreeParallelCellFinder
        (
            const fvPatch &
        );

        //- Construct from type and patch
        TreeParallelCellFinder
        (
            const word & TreeParallelCellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        TreeParallelCellFinder(const TreeParallelCellFinder &) = default;

        //- Destructor
        virtual ~TreeParallelCellFinder(){}

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
