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
    Foam::LawOfTheWall

Description
    Base class for laws of the wall.

Authors
    Timofey Mukha, Saleh Rezaeiravesh.

SourceFiles
    LawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef LawOfTheWall_H
#define LawOfTheWall_H

//#include "scalar.H"
//#include "dictionary.H"
//#include "word.H"
//#include "typeInfo.H"
//#include "runTimeSelectionTables.H"
//#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class LawOfTheWall
{
protected:
    
    //- Dictionary holding the model constants that the law uses
    dictionary constDict_;

public:

    TypeName("LawOfTheWall");
    
    // Constructors
    
        LawOfTheWall(const dictionary & dict)
        :
        constDict_(dict)
        {}
        
        LawOfTheWall(const word& lawName, const dictionary & dict)
        :
        constDict_(dict)
        {}

    // Destructor
        virtual ~LawOfTheWall() {}

    
    // RTS tables
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

        // RTS table "TypeAndDictionary"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            LawOfTheWall,
            TypeAndDictionary,
            (
                const word & lawName,
                const dictionary & dict
            ),
            (lawName, dict)
        )  
    
    // Selectors
        static autoPtr<LawOfTheWall> New 
        (
            const dictionary & dict
        );

        static autoPtr<LawOfTheWall> New 
        (
            const word & lawName,
            const dictionary & dict
        );
        
    // Member Functions
    
        //- Print info
        virtual void printCoeffs() const = 0; 

        //- Return the value of the implicit function defining the law
        virtual scalar
        value(scalar u, scalar y, scalar uTau, scalar nu) const = 0;

        //- Return the value of the  derivative of the implicit function
        //  defining the law
        virtual scalar 
        derivative(scalar u, scalar y, scalar uTau, scalar nu) const = 0;
        
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
