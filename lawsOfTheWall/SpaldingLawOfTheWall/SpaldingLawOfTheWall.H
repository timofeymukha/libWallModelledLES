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
    SpaldingLawOfTheWall

Description
    The law of the wall proposed by Spalding.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    SpaldingLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef SpaldingLawOfTheWall_H
#define SpaldingLawOfTheWall_H

#include "scalar.H"
#include "LawOfTheWall.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SpaldingLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class SpaldingLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa;
        
        //- The B model constant
        scalar B;
        
        //- The value of viscosity
        scalar nu;

public:

    // Constructors
        //- Construct given viscosity
        SpaldingLawOfTheWall(scalar nu);
    
        //- Construct provided kappa, B and viscosity
        SpaldingLawOfTheWall(scalar kappa, scalar B, scalar nu);

    // Destructor
        virtual ~SpaldingLawOfTheWall() {};

    // Member Functions

        //- Return the value of the implicit function defining the law
        virtual scalar value(scalar u, scalar y, scalar uTau) const;
        
        //- Return the value of the derivative of the implicit function
        //  defining the law
        scalar derivativeValue(scalar u, scalar y, scalar uTau) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
