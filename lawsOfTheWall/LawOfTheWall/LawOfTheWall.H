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
    LawOfTheWall

Description
    Base class for laws of the wall.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    LawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef LawOfTheWall_H
#define LawOfTheWall_H

#include "scalar.H"
#include "dictionary.H"
#include "word.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class LawOfTheWall
{

public:

    TypeName("LawOfTheWall");
    
    // Constructors
        LawOfTheWall() {}
    
        LawOfTheWall(const dictionary &) {}

    // Destructor
        virtual ~LawOfTheWall() {}

        
    // RTS table "Dictionary"    
    declareRunTimeSelectionTable
    (
        autoPtr,
        LawOfTheWall,
        Dictionary,
        (
            const dictionary & dict
        ),
        (dict)
    )    
    
    // Selectors
    static autoPtr<LawOfTheWall> New 
    (
        const dictionary & dict
    );
        
    // Member Functions
    
        //- Print the model coefficients
        virtual void printCoeffs() const = 0; 

        //- Return the value of the implicit function defining the law
        virtual scalar value(scalar u, scalar y, scalar uTau, scalar nu) const = 0;

        virtual scalar derivative(scalar u, scalar y, scalar uTau, scalar nu) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
