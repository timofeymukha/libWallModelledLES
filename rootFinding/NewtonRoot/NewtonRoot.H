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
    Newton root finder.
Author
    Timofey Mukha. All rights reserved.

SourceFiles
    NewtonRoot.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRoot_H
#define NewtonRoot_H

#include "scalar.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NewtonRoot Declaration
\*---------------------------------------------------------------------------*/

class NewtonRoot
{
    // Private data

        //- Pointer to the function
        std::function<scalar(scalar)> f_;

        //- Pointer to the derivative
        std::function<scalar(scalar)> d_;

        //- Tolerance
        const scalar eps_;
        
        //- Maximum number of iterations
        const label maxIter_;


public:

    // Static data members

        //- Maximum number of iterations
        static const label maxIter;


    // Constructors

        //- Construct given a function
        NewtonRoot(std::function<scalar(scalar)>,
                   std::function<scalar(scalar)>,
                   const scalar eps,
                   const label maxIter);


    // Destructor - default


    // Member Functions

        //- Return root
        scalar root(scalar guess) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
