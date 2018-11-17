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
    Foam::SampledField

Description
    Base class for a sampled field. Handles sampling of recomputing a 
    specific field. Also registers the field to the object registry
    under the sub-registry wallModelSampling -- patch-name.
    

Contributors/Copyright:
    2018 Timofey Mukha

SourceFiles
    SampledField.C

\*---------------------------------------------------------------------------*/

#ifndef SampledField_H
#define SampledField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                        Class SampledField Declaration
\*---------------------------------------------------------------------------*/

class SampledField
{
protected:
    
    // Protected data
        const fvPatch & patch_;
    
        const fvMesh & db_;
        
        const labelList cellIndexList_;
                    
public:

    // Constructors

        //- Construct from patch and cell index list
        SampledField
        (
            const fvPatch & patch,
            const labelList &  cil
        )
        :
            patch_(patch),
            db_(patch_.boundaryMesh().mesh()),
            cellIndexList_(cil)
        {}
        
        //- Clone the object
        virtual SampledField* clone() const = 0;

        //- Destructor
        virtual ~SampledField() {};        
            
    // Member functions
        
        //- Get the mesh
        const fvMesh & db() const
        {
            return db_;
        }
        
        //- Get the cell index list
        const labelList & cellIndexList() const
        {
            return cellIndexList_;
        }
    
        //- Get the patch
        const fvPatch & patch() const
        {
            return patch_;
        }
        
        //- Get the name of the sampled field
        virtual word name() const = 0;
        
        //- The number of dimensions of the field
        virtual label nDims() const = 0;
        
        //- Sample the field
        virtual void sample(scalarListList &) const = 0;
        
        //- Project vector field to wall-parallel direction
        void projectVectors(scalarListList &) const;

        void projectVectors(vectorField &) const;
        
        //- Register appropriate fields in the object registry
        virtual void registerFields() const {}
        
        //- Recompute the field
        virtual void recompute() const = 0;       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

