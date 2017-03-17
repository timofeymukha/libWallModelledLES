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
    WernerWengleLawOfTheWall

Description
    The law of the wall proposed by Werner and Wengle.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    WernerWengleLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef WernerWengleLawOfTheWall_H
#define WernerWengleLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class WernerWengleLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class WernerWengleLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar A_;
        
        //- The B model constant
        scalar B_;
        
public:
    
    TypeName("WernerWengle");

    // Constructors
        //- Construct with default values
    //    WernerWengleLawOfTheWall();
    
        //- Construct provided kappa, B and viscosity
      //  WernerWengleLawOfTheWall(scalar A, scalar B);
        
        //- Construct provided dictionary
        WernerWengleLawOfTheWall(const dictionary &);
        
        WernerWengleLawOfTheWall(const word &, const dictionary &);

    // Destructor
        virtual ~WernerWengleLawOfTheWall() {};

    // Member Functions
        
        scalar A() const
        {
            return A_;
        }
        
        scalar B() const
        {
            return B_;
        }
        
        //- Print info
        virtual void printCoeffs() const;

        //- Return the value of the implicit function defining the law
        virtual scalar value(scalar u, scalar y, scalar uTau, scalar nu) const;
        
        //- Return the value of the derivative of the implicit function
        //  defining the law
        virtual scalar
        derivative(scalar u, scalar y, scalar uTau, scalar nu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
