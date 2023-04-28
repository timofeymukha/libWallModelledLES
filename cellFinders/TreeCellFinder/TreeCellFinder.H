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
    Foam::TreeCellFinder

@brief
    Class for searching for sampling cells using.

Contributors/Copyright:
    2019-2021 Timofey Mukha


SourceFiles
    TreeCellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef TreeCellFinder_H
#define TreeCellFinder_H

#include "tmp.H"
#include "fixedValueFvPatchFields.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "CellFinder.H"
#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TreeCellFinder Declaration
\*---------------------------------------------------------------------------*/

class TreeCellFinder : public CellFinder
{

protected:
    
    // Protected data
        
    // Protected Member Functions

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("TreeCellFinder");
#endif

    // Constructors

        //- Construct from patch
        TreeCellFinder
        (
            const fvPatch &
        );

        //- Construct from type, patch
        TreeCellFinder
        (
            const word & TreeCellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        TreeCellFinder(const TreeCellFinder &) = default;
        
        //- Destructor
        virtual ~TreeCellFinder(){}

    // Member functions

        //- Compute distance field
        tmp<volScalarField> distanceField() const;

        //- Find the sampling cell indices for a single cell sampler
        void findCellIndices
        (
            labelList & indexList,
            const scalarField & h
        ) const;

        //- Find the sampling cell indices for a multi cell sampler
        void findCellIndices
        (
            labelListList & indexList,
            const scalarField & h,
            const bool excludeWallAdjacent
        ) const;
        
        //- Find cells close than 2max(h) to the wall
        tmp<labelField> findCandidateCellLabels
        (
            const scalarField & dist,
            const scalarField & h
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

