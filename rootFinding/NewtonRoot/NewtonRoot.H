/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    NewtontRoot

Description

Author
    Timofey Mukha. All rights reserved.

SourceFiles
    NewtonRoot.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRoot_H
#define NewtonRoot_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NewtonRoot Declaration
\*---------------------------------------------------------------------------*/

//template<typename Func>
class NewtonRoot
{
    // Private data

        //- Pointer to the function
        scalar (*f_) (scalar);

        //- Pointer to the derivative
        scalar (*d_) (scalar) ;

        //- Tolerance
        const scalar eps_;


public:

    // Static data members

        //- Maximum number of iterations
        static const label maxIter;


    // Constructors

        //- Construct given a function
        NewtonRoot(scalar (*f) (scalar), scalar (*d) (scalar), const scalar eps);


    // Destructor - default


    // Member Functions

        //- Return root
        scalar root(scalar xN) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NewtonRoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
