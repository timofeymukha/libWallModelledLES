/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
wallModel

Description
    Base abstract class for LES wall models.

Authors
    Timofey Mukha.  All rights reserved.

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
    
        //- The dictionary defining the boundary condition
        dictionary dict_;

    // Protected Member Functions

        //- Check the type of the patch
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
        //  nutWallFunctionFvPatchScalarField
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

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

