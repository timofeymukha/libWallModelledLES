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
    Foam::CellFinder

@brief
    Class for searching for sampling cells.

Contributors/Copyright:
    2019 Timofey Mukha


SourceFiles
    CellFinder.C

\*---------------------------------------------------------------------------*/

#ifndef CellFinder_H
#define CellFinder_H

#include "tmp.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CellFinder Declaration
\*---------------------------------------------------------------------------*/

class CellFinder
{

protected:
    
    // Protected data
        
        //- The patch to build the list for
        const fvPatch & patch_;

        //- The global object registry
        const fvMesh & mesh_;

    // Protected Member Functions

        void createFields() const;
        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("CellFinder");
#endif

    // Constructors

        //- Construct from patch and averaging time
        CellFinder
        (
            const fvPatch &
        );

        //- Construct from type, patch and averaging time
        CellFinder
        (
            const word & cellFinderName,
            const fvPatch &
        );

        //- Copy constructor
        CellFinder(const CellFinder &);
        
        //- Destructor
        virtual ~CellFinder();

    // Selectors
        static autoPtr<CellFinder> New
        (
            const word & CellFinderName,
            const fvPatch &
        );

        static autoPtr<CellFinder> New
        (
            const dictionary &,
            const fvPatch &
        );

    // Member functions

        //- Return the patch
        const fvPatch & patch() const
        {
            return patch_;
        }
        
        //- Get the mesh
        const fvMesh & mesh() const
        {
            return mesh_;
        }

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    // RTS tables
        
        // RTS table "Patch"
        declareRunTimeSelectionTable
        (
            autoPtr,
            CellFinder,
            Patch,
            (
                const word& CellFinderName,
                const fvPatch& p
            ),
            (CellFinderName, p)
        );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

