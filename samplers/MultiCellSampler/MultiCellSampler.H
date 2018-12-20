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
    Foam::MultiCellSampler

Description
    Class for sampling from several consequtive cells.

Contributors/Copyright:
    2018 Timofey Mukha

SourceFiles
    MultiCellSampler.C

\*---------------------------------------------------------------------------*/

#ifndef MultiCellSampler_H
#define MultiCellSampler_H

#include "fixedValueFvPatchFields.H"
#include "Sampler.H"
#include "SampledField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    
/*---------------------------------------------------------------------------*\
                         Class MultiCellSampler Declaration
\*---------------------------------------------------------------------------*/

class MultiCellSampler: public Sampler
{
protected:
    
    // Protected data
        
        //- The lists of indices of the cells that data is sampled from
        labelList indexList_;
        
        //- A list of wall-normal length-scales associated with the cells
        scalarField lengthList_;
        
        //- The distance from the wall that data is sampled from
        scalarField h_;
        
    // Protected Member Functions
                    
        //- Create list of cell-indices from where data is sampled
        virtual void createIndexList();
        
        //- Compute the length-scales
        virtual void createLengthList();
        
public:

    //- Runtime type information
        TypeName("MultiCellSampler");

    // Constructors

        //- Construct from patch and averaging time
        MultiCellSampler
        (
            const fvPatch&,
            scalar averagingTime
        );

        //- Construct from type, patch and averaging time
        MultiCellSampler
        (
            const word & samplerName,
            const fvPatch & p,
            scalar averagingTime
        );
        
        //- Copy constructor
        MultiCellSampler(const MultiCellSampler &) = default;
        
    // Destructor
        virtual ~MultiCellSampler();

    // Member functions

        //- Return the list of cell-indices that are used to sample data
        labelList indexList() const
        {
            return indexList_;
        }
        
        //- Return h
        const scalarField & h() const
        {
            return h_;
        }
                
        //- Return the length-list
        virtual const scalarField & lengthList() const
        {
            return lengthList_;
        }
        
        //- Element access operator
        inline label operator[](const label i) const
        {
            return indexList_[i];
        }
        
        //- Sample the fields
        virtual void sample() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

