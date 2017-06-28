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
    JohnsonKingEddyViscosity
Description
    The eddy viscosity model proposed by Johnson and King.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    JohnsonKingEddyViscosity.C
lablabla
\*---------------------------------------------------------------------------*/

#ifndef JohnsonKingEddyViscosity_H
#define JohnsonKingEddyViscosity_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "EddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class JohnsonKingEddyViscosity Declaration
\*---------------------------------------------------------------------------*/

class JohnsonKingEddyViscosity: public EddyViscosity
{
    
    // Private data
        
        //- The A+ model constant
        scalar APlus_;
        
        //- The kappa model constant
        scalar kappa_;
        
public:
    
    TypeName("JohnsonKing");

    // Constructors
        //- Construct with default values
    //    JohnsonKingEddyViscosity();
    
        //- Construct provided kappa, B and viscosity
      //  JohnsonKingEddyViscosity(scalar A, scalar B);
        
        //- Construct provided dictionary
        JohnsonKingEddyViscosity(const dictionary &);
        
        JohnsonKingEddyViscosity(const word &, const dictionary &);

    // Destructor
        virtual ~JohnsonKingEddyViscosity() {};

    // Member Functions
        
        scalar APlus() const
        {
            return APlus_;
        }
        
        scalar kappa() const
        {
            return kappa_;
        }
        
        //- Print info
        virtual void printCoeffs() const;

        //- Return the value of nut
        virtual scalarList
        value(const scalarList & y, scalar uTau, scalar nu) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
