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
    Foam::wallModelFvPatchScalarField

Description
    Base abstract class for LES wall models.
    
    Handles creating and storing fields used by the wall models.
    data. The following fields are created and stored in the registry.
    - h, the patch fields of which hold the distance to the cells used for
    sampling. This field is read and must be present for the simulation to run.
    - samplingCells, which marks the cells used for sampling with a value
    corresponding to the index of the patch that uses them.
    - uTauPredicted, the patch fields of which hold the value of uTau as
    predicted by the wall model.
    - wallShearStress, the patch fields of which hold the value of the wall 
    shear stress.
    - wallGradU, the patch fields of which store the wall-nromal of the velocity
    gradient.


Contributors/Copyright:
    2018-2019 Timofey Mukha

SourceFiles
    wallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallModelFvPatchScalarField_H
#define wallModelFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Sampler.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class wallModelFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallModelFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    //- Fraction of CPU time spent on wall modelling
    scalar consumedTime_;
    
protected:

    // Protected data
           
        //- Timescale of the averaging
        scalar averagingTime_;

        //- Create fields and add to registry
        virtual void createFields() const;

    // Protected Member Functions

        //- Check that the patch is a wall
        virtual void checkType();

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const = 0;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("wallModel");

    // Constructors

        //- Construct from patch and internal field
        wallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wallModelFvPatchScalarField
        //  onto a new patch
        wallModelFvPatchScalarField
        (
            const wallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallModelFvPatchScalarField
        (
            const wallModelFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        wallModelFvPatchScalarField
        (
            const wallModelFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions
      
        scalar averagingTime() const
        {
            return averagingTime_;
        }

        scalar consumedTime() const
        {
            return consumedTime_;
        }
        
        //- Update the boundary values
        virtual void updateCoeffs();

        //- Write to stream
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

