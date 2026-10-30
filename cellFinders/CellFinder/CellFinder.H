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

    Derived classes are selected at run time by the value of the sampler
    keyword, with CellFinder appended, e.g. Crawling selects the
    CrawlingCellFinder.

    The common interface, findCells(), returns the sampling locations found
    on this processor together with a map, which distributes values computed
    at these locations to the faces of the patch. Cell finders that only
    search on the processor holding the wall face return a map without
    communication, while others may return sampling locations for the wall
    faces of other processors.

Contributors/Copyright:
    2019-2026 Timofey Mukha


SourceFiles
    CellFinder.C
    CellFinderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef CellFinder_H
#define CellFinder_H

#include "tmp.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CellFinder Declaration
\*---------------------------------------------------------------------------*/

class CellFinder
{

public:

    //- A sampling location for a wall face. Also used to pass the search
    //- state between the processors.
    class samplingLocation
    {
    public:

        //- The processor holding the wall face
        label originProc;

        //- The index of the wall face in the patch on originProc
        label originFace;

        //- The cell, local to the current processor
        label cell;

        //- The face through which the cell is entered when crawling. When
        //- passed between processors, the index of the face in the
        //- processor patch.
        label face;

        //- The consecutive off-wall index of the cell, starting from 1
        label layer;

        //- The centre of the wall face
        point wallCentre;

        //- The outward unit normal of the wall face
        vector wallNormal;

        //- The prescribed sampling height
        scalar h;

        friend bool operator==
        (
            const samplingLocation & a,
            const samplingLocation & b
        )
        {
            return
                a.originProc == b.originProc
             && a.originFace == b.originFace
             && a.cell == b.cell
             && a.face == b.face
             && a.layer == b.layer
             && a.wallCentre == b.wallCentre
             && a.wallNormal == b.wallNormal
             && a.h == b.h;
        }

        friend bool operator!=
        (
            const samplingLocation & a,
            const samplingLocation & b
        )
        {
            return !(a == b);
        }

        friend Ostream & operator<<(Ostream & os, const samplingLocation & l)
        {
            os  << l.originProc << token::SPACE << l.originFace
                << token::SPACE << l.cell << token::SPACE << l.face
                << token::SPACE << l.layer << token::SPACE << l.wallCentre
                << token::SPACE << l.wallNormal << token::SPACE << l.h;
            return os;
        }

        friend Istream & operator>>(Istream & is, samplingLocation & l)
        {
            is  >> l.originProc >> l.originFace >> l.cell >> l.face
                >> l.layer >> l.wallCentre >> l.wallNormal >> l.h;
            return is;
        }
    };

protected:
    
    // Protected data
        
        //- The patch to build the list for
        const fvPatch & patch_;

        //- The global object registry
        const fvMesh & mesh_;

    // Protected Member Functions

        void createFields() const;

        //- Create a sampling location in the wall-adjacent cell of a face
        //- of the patch
        samplingLocation wallLocation
        (
            const label patchFaceI,
            const scalar h,
            const vectorField & faceNormals
        ) const;

        //- Create sampling locations on this processor for the given cells,
        //- one per face of the patch, and the corresponding map, which does
        //- not involve communication
        void createLocalLocations
        (
            const labelList & cells,
            const scalarField & h,
            List<samplingLocation> & locations,
            mapDistributeBase & map
        ) const;

        //- Construct the map from sampling locations, which may correspond
        //- to faces on any processor, to the faces of the patch. Must be
        //- called on all processors.
        void createMap
        (
            const UList<samplingLocation> & locations,
            mapDistributeBase & map
        ) const;

        //- Send each list to the respective processor and receive the lists
        //- sent to this processor. Must be called on all processors.
        template<class Type>
        static void exchange
        (
            const UList<DynamicList<Type>> & sendLists,
            List<List<Type>> & recvLists
        );

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

        //- Find a sampling cell for each face of the patch. On output, the
        //- sampling locations on this processor and the map, which
        //- distributes values computed at them to the faces of the patch.
        //- Must be called on all processors.
        virtual void findCells
        (
            List<samplingLocation> & locations,
            mapDistributeBase & map,
            const scalarField & h,
            const bool hIsIndex
        ) const = 0;

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

#ifdef NoRepository
    #include "CellFinderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

