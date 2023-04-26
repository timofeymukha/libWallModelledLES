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
    Foam::EddyViscosity

@brief
    Base abstract class for 1D turbulence models based on eddy viscosity.

    Meant to be used with ODE based wall models that need the distribution
    of eddy viscosity on a wall-normal 1D mesh associated with each face.

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    LawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef EddyViscosity_H
#define EddyViscosity_H

#include "refCount.H"
#include "Field.H"
#include "SingleCellSampler.H"

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

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("EddyViscosity");
#endif
    
    // Constructors
    
        EddyViscosity(const dictionary & dict)
        :
        constDict_(dict)
        {}
        
        EddyViscosity
        (
            const word& modelName,
            const dictionary & dict
        )
        :
        constDict_(dict)
        {}

        EddyViscosity() = default;

        //- Copy constructor
        EddyViscosity(const EddyViscosity & orig) = default;

        //- Assignemnt
        EddyViscosity & operator=(const EddyViscosity &) = default;

        //- Clone
        virtual autoPtr<EddyViscosity> clone() const = 0;

    // Destructor
        virtual ~EddyViscosity() {}

    
    // Selectors
        static autoPtr<EddyViscosity> New 
        (
            const dictionary & dict
        );

        static autoPtr<EddyViscosity> New 
        (
            const word & modelName,
            const dictionary & dict
        );
        
    // Member Functions

        //- Add necessary sampled fields to the sampler
        virtual void addFieldsToSampler(Sampler & sampler) const
        {}
    
        //- Print info
        virtual void printCoeffs() const = 0; 

        //- Return the value of nut
        virtual scalarList
        value
        (
            const SingleCellSampler & sampler,
            const label index,
            const scalarList & y,
            const scalar uTau,
            const scalar nu
        ) const = 0;
        
        //- Write information about the law to stream
        virtual void write(Ostream & os) const; 
        
        //- Get the dictionary with the constants
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
            EddyViscosity,
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
            EddyViscosity,
            TypeAndDictionary,
            (
                const word & modelName,
                const dictionary & dict
            ),
            (modelName, dict)
        ); 
#endif

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
