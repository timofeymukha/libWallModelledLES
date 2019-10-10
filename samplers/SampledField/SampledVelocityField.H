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
    Foam::SampledVelocityField

Description
    Class for the sampled velocity field.

Contributors/Copyright:
    2019 Timofey Mukha

SourceFiles
    SampledVelocityField.C

\*---------------------------------------------------------------------------*/

#ifndef SampledVelocityField_H
#define SampledVelocityField_H

#include "fixedValueFvPatchFields.H"
#include "SampledField.H"
#include "scalarListListIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SampledVelocityField Declaration
\*---------------------------------------------------------------------------*/

class SampledVelocityField
:
    public SampledField
{

public:

    // Constructors

        //- Construct from patch and interpolation type
        SampledVelocityField
        (
            const fvPatch & patch,
            word interpolationType = "cell"
        )
        :
            SampledField(patch, interpolationType)
        {
            createField();
            setInterpolator(interpolationType);
        }

        SampledVelocityField(const SampledVelocityField & orig)
        :
            SampledField(orig)
        {
            if (orig.interpolator_.get() != nullptr)
            {
                setInterpolator(orig.interpolator_().type())
            }
        }
        //- Clone the object
        autoPtr<SampledField> clone() const override
        {
            return autoPtr<SampledField>
            (
                new SampledVelocityField(*this)
            );
        }
        
        //- Destructor
        virtual ~SampledVelocityField() {}
            
    // Member functions
        
        //- Sample the velocity values
        void sample
        (
            scalarListList & sampledValues,
            const labelList & indexList,
            const scalarField & h
        ) const override;

        //- Sample the velocity value from multile cellss
        void sample
        (
            scalarListListList &,
            const labelListList &
        ) const override;
                
        //- The number of dimensions of the field
        label nDims() const override
        {
            return 3;
        }
        
        //- Get the name of the sampled field
        word name() const override
        {
            return "U";
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
        void recompute() const override {}

        //- Set the inerpolator
        void setInterpolator(const word interpolationType) override;

        //- Create the field. Redunant, U should be created by the solver.
        void createField() const override {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

