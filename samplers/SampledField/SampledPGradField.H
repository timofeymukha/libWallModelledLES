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
    Foam::SampledPGradField

Description
    Class for a sampled pressure gradient field.

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    SampledPGradField.C

\*---------------------------------------------------------------------------*/

#ifndef SampledPGradField_H
#define SampledPGradField_H

#include "SampledField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                         Class SampledPGradField Declaration
\*---------------------------------------------------------------------------*/

class SampledPGradField
:
    public SampledField
{
              
public:

    // Constructors

        //- Construct from patch and averaging time
        
        SampledPGradField
        (
            const fvMesh & db
        )
        :
            SampledField(db)
        {
            registerFields();
        }

        virtual ~SampledPGradField() {};        
            
    // Member functions

        // Register appropriate fields in the object registry
        virtual void registerFields() const;
        
        //- Recompute field
        virtual void recompute() const;       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

