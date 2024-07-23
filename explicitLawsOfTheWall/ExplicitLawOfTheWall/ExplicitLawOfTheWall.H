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
    Foam::ExplicitLawOfTheWall

@brief
    Base abstract class for explicit laws of the wall.

    It is to be used in conjuction with the ExplicitLOTW wall model.
     The law of the wall provides an expression

Contributors/Copyright:
    2024 Timofey Mukha


\*---------------------------------------------------------------------------*/

#ifndef ExplicitLawOfTheWall_H
#define ExplicitLawOfTheWall_H

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
                      Class ExplicitLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class ExplicitLawOfTheWall : public refCount
{
protected:

    //- Dictionary holding the model constants that the law uses
    dictionary constDict_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("ExplicitLawOfTheWall");
#endif

    // Constructors

        //- Construct from dictionary
        ExplicitLawOfTheWall(const dictionary & dict)
        :
        constDict_(dict)
        {}

        //- Construct from TypeName and dictionary
        ExplicitLawOfTheWall
        (
            const word& lawName,
            const dictionary & dict
        )
        :
        constDict_(dict)
        {}

        //- Default constructor
        ExplicitLawOfTheWall() = default;

        //- Assignment
        ExplicitLawOfTheWall & operator=(const ExplicitLawOfTheWall &) = default;

        //- Copy constructor
        ExplicitLawOfTheWall(const ExplicitLawOfTheWall & orig)
        :
        refCount(),
        constDict_(orig.constDict_)
        {}

        //- Destructor
        virtual ~ExplicitLawOfTheWall() {}

        //- Clone
        virtual autoPtr<ExplicitLawOfTheWall> clone() const = 0;

    // Selectors
        static autoPtr<ExplicitLawOfTheWall> New
        (
            const dictionary & dict
        );

        static autoPtr<ExplicitLawOfTheWall> New
        (
            const word & lawName,
            const dictionary & dict
        );

    // Member Functions

        //- Add necessary sampled fields to the sampler
        virtual void addFieldsToSampler(Sampler & sampler){}

        //- Print info to terminal
        virtual void printCoeffs() const = 0;

        //- Return the value of the friction velocity
        virtual scalar uTau
        (
            const SingleCellSampler & sampler,
            label index,
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
            ExplicitLawOfTheWall,
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
            ExplicitLawOfTheWall,
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
