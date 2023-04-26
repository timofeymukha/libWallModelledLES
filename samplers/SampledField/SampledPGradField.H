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

@brief
    Class for the sampled pressure gradient field.

Contributors/Copyright:
    2018-2023 Timofey Mukha

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
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("SampledPGradField");
#endif

    // Constructors

        //- Construct from patch and interpolation type
        SampledPGradField
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
        SampledPGradField(const SampledPGradField & orig)
        :
            SampledField(orig)
        {}

        //- Clone the object
        autoPtr<SampledField> clone() const override
        {
            return autoPtr<SampledField>
            (
                new SampledPGradField(*this)
            );
        }
        
        //- Destructor
        virtual ~SampledPGradField() {}
            
    // Member functions
        
        //- The number of dimensions of the field
        label nDims() const override
        {
            return 3;
        }
        
        //- Get the name of the sampled field
        word name() const override
        {
            return "pGrad";
        }

        //- Sample the pressure gradient values from a single cell per wall face
        void sample
        (
            scalarListList & sampledValues,
            const labelList & indexList,
            const scalarField & h
        ) const override;

        //- Sample the pressure gradient values from multiple cells per wall face
        void sample(scalarListListList &, const labelListList &) const override;
        
        //- Register the pGrad field in the object registry
        void registerFields
        (
            const labelList &
        ) const override;

        //- Register the pGrad field in the object registry
        void registerFields
        (
            const labelListList &
        ) const override;

        //- Recompute the field
        void recompute() const override;

        //- Create the global pressure gradient field that will be sampled
        virtual void createField() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

