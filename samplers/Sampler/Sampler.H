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

Contributors/Copyright:
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh


SourceFiles
    Sampler.C

\*---------------------------------------------------------------------------*/

#ifndef Sampler_H
#define Sampler_H

#include "fixedValueFvPatchFields.H"
#include "SampledField.H"


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
        
        //- The global object registry
        const fvMesh & mesh_;
        
        //- Whether we need to compute the pressure gradient
        List<SampledField *> sampledFields_;
                      
    // Protected Member Functions
                    
        //- Create list of cell-indices from where data is sampled
        void createIndexList();
        
        //- Compute the length-scales
        void createLengthList();
        
        //- Project a vector field
        void project(vectorField & field) const;
        
public:

    //- Runtime type information
    TypeName("Sampler");

    // Constructors

        //- Construct from patch and averaging time
        Sampler
        (
            const fvPatch&,
            scalar averagingTime
        );
        
        Sampler(const Sampler &);
        
    // Destructor
        virtual ~Sampler();

    // Member functions

        //- Create fields
        void createFields();

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
            return mesh().subRegistry("wallModelSampling").subRegistry(patch().name());
        }
        
        //- Get the mesh
        const fvMesh & mesh() const
        {
            return mesh_;
        }
        
        //- Element access operator
        inline label operator[](const label i) const
        {
            return indexList_[i];
        }
        
        //- Recompute fields to be sampled
        void recomputeFields() const;
        
        //- Sample the fields
        void sample() const; 
        
        //- Add field for sampling
        void addField(SampledField *);
        
        //- Convert sampled listList to a field
        template<class Type>
        void listListToField(const scalarListList &, Field<Type> &) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

