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
    intended for meshes with reasonably aligned near-wall layers, where
    following opposite faces approximates the wall-normal direction.

    The sampling height \f$h\f$ can be interpreted either as a physical
    distance from the wall face or as a consecutive off-wall cell index. In
    distance mode, crawling stops when the next opposing face is farther from
    the wall face than \f$h\f$, and the current cell is used. In index mode,
    \f$h\f$ is rounded to the requested layer index.

    For multi-cell sampling, every crawled cell is retained up to the requested
    distance or layer index. If the crawl reaches a boundary patch before the
    requested target is reached, the last valid cells are retained. The optional
    \c excludeWallAdjacent argument removes the wall-adjacent cell only when at
    least one further cell remains.

    Invalid distances or indices fall back to the wall-adjacent cell.

Contributors/Copyright:
    2019-2026 Timofey Mukha


SourceFiles
    CrawlingCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef CrawlingCellFinder_H
#define CrawlingCellFinder_H

#include "tmp.H"
#include "fixedValueFvPatchFields.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
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
    
    // Protected data
        
    // Protected Member Functions

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
        CrawlingCellFinder(const CrawlingCellFinder &) = default;
        
        //- Destructor
        virtual ~CrawlingCellFinder(){}

    // Member functions

        void findCellIndices
        (
            labelList & indexList,
            const scalarField & h,
            const bool hIsIndex
        ) const;

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
