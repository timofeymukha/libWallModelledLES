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

@brief
    Base class for a sampled field. Handles sampling of and recomputing a 
    specific field. Also registers the field to the object registry
    under the sub-registry wallModelSampling -- patch-name.
    

Contributors/Copyright:
    2018-2021 Timofey Mukha

SourceFiles
    SampledField.C

\*---------------------------------------------------------------------------*/

#ifndef SampledField_H
#define SampledField_H

#include "fixedValueFvPatchFields.H"
#include "scalarListListIOList.H"
#include "runTimeSelectionTables.H"
#include "interpolation.H"

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
        //- The patch
        const fvPatch & patch_;
    
        //- The mesh
        const fvMesh & mesh_;

        //- Interpolation type
        const word interpolationType_;
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("SampledField");
#endif

    // Constructors

        //- Construct from patch and interpolation type
        SampledField
        (
            const fvPatch & patch,
            const word interpolationType="cell"
        )
        :
            patch_(patch),
            mesh_(patch_.boundaryMesh().mesh()),
            interpolationType_(interpolationType)
        {
        }
      
        //- Copy constructor
        SampledField
        (
            const SampledField & orig
        )
        :
            patch_(orig.patch()),
            mesh_(orig.mesh()),
// #ifdef FOAM_AUTOPTR_HAS_CLONE_METHOD
            interpolationType_(orig.interpolationType())
// #else
//             interpolator_(orig.interpolator_, false)
// #endif
        {}

        //- Clone the object
        virtual autoPtr<SampledField> clone() const = 0;

        //- Destructor
        virtual ~SampledField() {}
    // Member functions
        
        //- Get the mesh
        const fvMesh & mesh() const
        {
            return mesh_;
        }

        //- Subregistry for sampled fields
        const objectRegistry & db() const
        {
            return mesh().subRegistry
            (
                "wallModelSampling"
            ).subRegistry(patch().name());
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
        virtual void sample
        (
            scalarListList & sampledValues,
            const labelList & indexList,
            const scalarField & h
        ) const = 0;

        //- Sample the field from multiple cells
        virtual void sample
        (
            scalarListListList & sampledValues,
            const labelListList & indexListList
        ) const = 0;
        
        //- Register appropriate fields in the object registry
        virtual void registerFields
        (
            const labelList &
        ) const = 0; 

        virtual void registerFields
        (
            const labelListList &
        ) const = 0; 
        
        //- Recompute the values of the field
        virtual void recompute() const = 0;

        //- Create the global field that will be sampled
        virtual void createField() const = 0;
        
        //- Get interpolation type
        word interpolationType() const
        {
            return interpolationType_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

