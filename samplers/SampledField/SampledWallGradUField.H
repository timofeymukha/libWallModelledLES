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
    Foam::SampledWallGradUField

@brief
    Class for the sampled velocity gradient field at the wall.

Contributors/Copyright:
    2019-2023 Timofey Mukha

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
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("SampledWallGradUField");
#endif

    // Constructors

        //- Construct from patch
        SampledWallGradUField
        (
            const fvPatch & patch,
            const word interpolationType="cell"
        )
        :
            SampledField(patch, interpolationType)
        {
            createField();
        }

        //- Copy constructor
        SampledWallGradUField(const SampledWallGradUField & orig)
        :
            SampledField(orig)
        {}

        //- Clone the object
        autoPtr<SampledField> clone() const override
        {
            return autoPtr<SampledField>
            (
                new SampledWallGradUField(*this)
            );
        }
        
        //- Destructor
        virtual ~SampledWallGradUField() {}
            
    // Member functions
        
        //- Sample the wall-normal velocity gradient values
        void sample
        (
            scalarListList & sampledValues,
            const labelList & indexList,
            const scalarField & h
        ) const override;

        //- Sample the wall-normal velocity gradient values from multiple cells
        void sample(scalarListListList &, const labelListList &) const override;
                
        //- The number of dimensions of the field
        label nDims() const override
        {
            return 3;
        }
        
        //- Get the name of the sampled field
        word name() const override
        {
            return "wallGradU";
        }
        
        //- Register appropriate fields in the object registry
        void registerFields
        (
            const labelList &
        ) const override;

        //- Register appropriate fields in the object registry
        void registerFields
        (
            const labelListList &
        ) const override;

        //- Recompute the field
        void recompute() const override;

        //- Create the global wallGrad gradient field that will be sampled
        virtual void createField() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

