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

@brief
    Base abstract class for laws of the wall.
 
    It is to be used in conjuction
    with the LOTW wall model. The law of the wall provides an implicit function
    \f$F(u, y, u_\tau, \nu)\f$ and its derivative, which can be used to
    iteratively solve for the friction velocity.

Authors
    Timofey Mukha, Saleh Rezaeiravesh.

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh


\*---------------------------------------------------------------------------*/

#ifndef LawOfTheWall_H
#define LawOfTheWall_H

#include "dictionary.H"
#include "refCount.H"
#include "word.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Sampler;
class SingleCellSampler;

/*---------------------------------------------------------------------------*\
                      Class LawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class LawOfTheWall : public refCount
{
protected:
    
    //- Dictionary holding the model constants that the law uses
    dictionary constDict_;
    
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("LawOfTheWall");
#endif
    
    // Constructors
    
        //- Construct from dictionary
        LawOfTheWall(const dictionary & dict)
        :
        constDict_(dict)
        {}
        
        //- Construct from TypeName and dictionary
        LawOfTheWall
        (
            const word& lawName,
            const dictionary & dict
        )
        :
        constDict_(dict)
        {}

        //- Default constructor
        LawOfTheWall() = default;

        //- Assignment
        LawOfTheWall & operator=(const LawOfTheWall &) = default;

        //- Copy constructor
        LawOfTheWall(const LawOfTheWall & orig)
        :
        refCount(),
        constDict_(orig.constDict_)
        {}

        //- Destructor
        virtual ~LawOfTheWall() {}

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const = 0;

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
    
        //- Add necessary sampled fields to the sampler
        virtual void addFieldsToSampler(Sampler & sampler){}

        //- Print info to terminal
        virtual void printCoeffs() const = 0; 

        //- Return the value of the implicit function defining the law
        virtual scalar value
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const = 0;

        //- Return the value of the  derivative of the implicit function
        //  defining the law. Should throw an error if derivative does not 
        //  exist.
        virtual scalar derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const = 0;
        
        //- Write information about the law to stream
        virtual void write(Ostream & os) const; 
        
        //- Return the dictionary with the constants
        dictionary constDict() const
        {
            return constDict_;
        }
        
        
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
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
        );

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
        );  
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
