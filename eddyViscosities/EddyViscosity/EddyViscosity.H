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
    EddyViscosity

Description
    Base class for 1d "mixing length"-type turbulence models.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    EddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef EddyViscosity_H
#define EddyViscosity_H

#include "scalar.H"
#include "dictionary.H"
#include "word.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "Field.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class EddyViscosity Declaration
\*---------------------------------------------------------------------------*/

class EddyViscosity
{
protected:
    
    //- Dictionary holding the model constants
    dictionary constDict_;

public:

    TypeName("EddyViscosity");
    
    // Constructors
    
        EddyViscosity(const dictionary & dict)
        :
        constDict_(dict)
        {}
        
        EddyViscosity(const word& modelName, const dictionary & dict)
        :
        constDict_(dict)
        {}

    // Destructor
        virtual ~EddyViscosity() {}

    
    // RTS tables
        // RTS table "Dictionary"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            EddyViscosity,
            Dictionary,
            (
                const dictionary & dict
            ),
            (dict)
        )

        // RTS table "TypeAndDictionary"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            EddyViscosity,
            TypeAndDictionary,
            (
                const word & modelName,
                const dictionary & dict
            ),
            (modelName, dict)
        )  
    
    // Selectors
        static autoPtr<EddyViscosity> New 
        (
            const dictionary & dict
        );

        static autoPtr<EddyViscosity> New 
        (
//            const word & lawName,
            const word & modelName,
            const dictionary & dict
        );
        
    // Member Functions
    
        //- Print info
        virtual void printCoeffs() const = 0; 

        //- Return the value of nut
        virtual scalarList
        value(const scalarList & y, scalar uTau, scalar nu) const = 0;
        
        //- Write information about the law to stream
        virtual void write(Ostream & os) const; 
        
        //- Get the dictionary with the constants
        dictionary constDict() const
        {
            return constDict_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
