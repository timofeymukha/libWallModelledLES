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
    Base class for a sampled field.

Authors
    Timofey Mukha, Saleh Rezaeiravesh

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
        
        const fvMesh & db_;
              
public:


    // Constructors

        //- Construct from patch and averaging time
        
        SampledField
        (
            const fvMesh & db
        )
        :
            db_(db)
        {}

        virtual ~SampledField() {};        
            
    // Member functions
        
        //- Get the mesh
        const fvMesh & db() const
        {
            return db_;
        }

        // Register appropriate fields in the object registry
        virtual void registerFields() const {}
        
        //- Recompute field
        virtual void recompute() const = 0;       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

