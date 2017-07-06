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
EquilibriumODEWallModel

Description
    ODE wall model with a zero source term.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    EquilibriumODEWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef EquilibriumODEWallModelFvPatchScalarField_H
#define EquilibriumODEWallModelFvPatchScalarField_H

#include "ODEWallModelFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class EquilibriumODEWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class EquilibriumODEWallModelFvPatchScalarField
:
    public ODEWallModelFvPatchScalarField
{
protected:
        
    // Protected Member Functions
        //- Write root finder and EquilibriumODE properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        
//        scalar source() const
//        scalar source(const int &) const 
        virtual void source(const int, const vector , vector &) const
        {};
        
public:

    //- Runtime type information
    TypeName("EquilibriumODEWallModel");


    // Constructors

        //- Construct from patch and internal field
        EquilibriumODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        EquilibriumODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  EquilibriumODEWallModelFvPatchScalarField
        //  onto a new patch
        EquilibriumODEWallModelFvPatchScalarField
        (
            const EquilibriumODEWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        EquilibriumODEWallModelFvPatchScalarField
        (
            const EquilibriumODEWallModelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new EquilibriumODEWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        EquilibriumODEWallModelFvPatchScalarField
        (
            const EquilibriumODEWallModelFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new EquilibriumODEWallModelFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
            
        
        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
