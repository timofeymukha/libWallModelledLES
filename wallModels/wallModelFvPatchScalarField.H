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
    
    Handles creating and storing fields
    used by the wall models as well as figuring what cells are used for sampling
    data. The following fields are created and stored in the registry.
    - h, the patch fields of which hold the distance to the cells used for
    sampling. This field is read and must be present for the simulation to run.
    - samplingCells, which marks the cells used for sampling with a value
    corresponding to the index of the patch that uses them.
    - uTau, the patch fields of which hold the value of uTau as predicted by the
    wall model.
    - magGradU, the patch fields of which store the magnitude of the velocity
    gradient.


Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    wallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallModelFvPatchScalarField_H
#define wallModelFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

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
protected:
    
    // Protected data
    
        //- The indices of the cells that data is sampled from
        labelList cellIndexList_;
        
        //- The distance from the wall that data is sampled from
        scalarField h_;
        
        //- The velocity values the wall model uses
        vectorField U_;
          
        //- Create fields and add to registry
        void createFields() const;

    // Protected Member Functions

        //- Check that the patch is a wall
        virtual void checkType();

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const = 0;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Benchmark for calculating the friction velocity
        //  using the default Spalding wall model
        tmp<scalarField> calcUTauBench(const scalarField& magGradU) const;

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

        //- Create list of cell-indices from where data is sampled
        void createCellIndexList();
        
        //- Sample velocity values from the LES
        void sample();
        
        //- Return the list of cell-indices that are used to sample data
        labelList cellIndexList() const
        {
            return cellIndexList_;
        }
        
        //- Return h
        scalarField h() const
        {
            return h_;
        }
        
        //- return U_
        vectorField U() const
        {
            return U_;
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

