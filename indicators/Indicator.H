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
    Foam::Indicator

@brief
    Base class for determining whether the TBL is turbulent.

Contributors/Copyright:
    2024 Timofey Mukha


SourceFiles
    Indicator.C

\*---------------------------------------------------------------------------*/

#ifndef Indicator_H
#define Indicator_H

#include "tmp.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Indicator Declaration
\*---------------------------------------------------------------------------*/

class Indicator
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
        TypeName("Indicator");
#endif
    // Constructors

        //- Construct from patch and averaging time
        Indicator
        (
            const fvPatch &
        );

        //- Construct from type, patch and averaging time
        Indicator
        (
            const word & indicatorName,
            const fvPatch &
        );

        //- Copy constructor
        Indicator(const Indicator &);

    // Selectors
        static autoPtr<Indicator> New
        (
            const word & IndicatorName,
            const fvPatch &
        );

        static autoPtr<Indicator> New
        (
            const dictionary &,
            const fvPatch &
        );

        //- Destructor
        virtual ~Indicator();

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

        //- Compute the indicator field
        void compute(const volScalarField & nu) const;

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    // RTS tables

        // RTS table "Patch"
        declareRunTimeSelectionTable
        (
            autoPtr,
            Indicator,
            Patch,
            (
                const word& IndicatorName,
                const fvPatch& p
            ),
            (IndicatorName, p)
        );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

