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
    Foam::Sampler

Description
    Class for sampling data to the wall models.

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    Sampler.C

\*---------------------------------------------------------------------------*/

#ifndef Sampler_H
#define Sampler_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Sampler Declaration
\*---------------------------------------------------------------------------*/

class Sampler
{
protected:
    
    // Protected data
        
        //- The patch to build the list for
        const fvPatch & patch_;    
    
        //- The indices of the cells that data is sampled from
        labelList indexList_;
        
        //- A list of wall-normal length-scales associated with the cells
        scalarField lengthList_;
        
        //- The distance from the wall that data is sampled from
        scalarField h_;
        
        //- Time-averaging scale of the sampled values
        scalar averagingTime_;
        
        //- The object registry for sampled data
        const objectRegistry & db_;
        
        //- The global object registry
        const objectRegistry & mesh_;
                
    // Protected Member Functions
        
        //- Create h field and add to registry
        void createFields() const;
             
        //- Create list of cell-indices from where data is sampled
        void createIndexList();
        
        //- Compute the length-scales
        void createLengthList();
        
        //- Project a vector field
        void project(vectorField & field) const;
        
public:


    // Constructors

        //- Construct from patch and averaging time
        Sampler
        (
            const fvPatch&,
            scalar averagingTime
        );

    // Member functions
      
        //- Return the list of cell-indices that are used to sample data
        labelList indexList() const
        {
            return indexList_;
        }
        
        //- Return h
        scalarField h() const
        {
            return h_;
        }
                
        //- Return the length-list
        scalarField lengthList() const
        {
            return lengthList_;
        }
        
        //- Return the patch
        const fvPatch & patch() const
        {
            return patch_;
        }
        
        //- Get the registry with sampled fields
        const objectRegistry & db() const
        {
            return db_;
        }
        
        //- Get the mesh
        const objectRegistry & mesh() const
        {
            return mesh_;
        }
        
        //- Element access operator
        inline label operator[](const label i) const
        {
            return indexList_[i];
        }
        
        //- Sample the fields
        void sample() const;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
