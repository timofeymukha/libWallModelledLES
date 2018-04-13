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
    Class for the sampled velocity field.

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    SampledWallGradUField.C

\*---------------------------------------------------------------------------*/

#ifndef SampledWallGradUField_H
#define SampledWallGradUField_H

#include "fixedValueFvPatchFields.H"
#include "SampledField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                    Class SampledWallGradUField Declaration
\*---------------------------------------------------------------------------*/

class SampledWallGradUField
:
    public SampledField
{
   
public:
    
    // Constructors

        //- Construct from patch and cell index list
        SampledWallGradUField
        (
            const fvPatch & patch,
            const labelList & cil    
        )
        :
            SampledField(patch, cil)
        {
            registerFields();
        }
            
        //- Clone the object    
        virtual SampledWallGradUField* clone() const
        {
            return new SampledWallGradUField(*this);
        }
        
        //- Destructor
        virtual ~SampledWallGradUField() {};        
            
    // Member functions
        
        //- Sample the wall-normal velocity gradient values
        virtual void sample(scalarListList &) const;
                
        //- The number of dimensions of the field
        virtual label nDims() const
        {
            return 3;
        }
        
        //- Get the name of the sampled field
        virtual word name() const
        {
            return "wallGradU";
        }
        
        //- Register appropriate fields in the object registry
        virtual void registerFields() const;
        
        //- Recompute the field
        virtual void recompute() const;   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

