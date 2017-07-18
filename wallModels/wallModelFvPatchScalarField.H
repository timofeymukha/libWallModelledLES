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
wallModel

Description
    Base abstract class for LES wall models.

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
          
        //- Create fields and add to registry
        void createFields() const;

    // Protected Member Functions

        //- Check the type of the patch
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
        
        // Evaluation functions

            //- Update the boundary values
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

