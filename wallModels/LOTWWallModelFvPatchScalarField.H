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
LOTWWallModel

Description
    Class for the wall models based on a Law of The Wall.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    LOTWWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef LOTWWallModelFvPatchScalarField_H
#define LOTWWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class LOTWWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class LOTWWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Data
        autoPtr<RootFinder> rootFinder_;
    
        //- The "law of the wall" to be used
        autoPtr<LawOfTheWall> law_;
    
    // Protected Member Functions
        //- Write root finder and law properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;
        
        //- Benchmark for calculating the friction velocity
        virtual tmp<scalarField> calcUTauBench(const scalarField& magGradU) const;
        



public:

    //- Runtime type information
    TypeName("LOTWWallModel");


    // Constructors

        //- Construct from patch and internal field
        LOTWWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        LOTWWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  LOTWWallModelFvPatchScalarField
        //  onto a new patch
        LOTWWallModelFvPatchScalarField
        (
            const LOTWWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        LOTWWallModelFvPatchScalarField
        (
            const LOTWWallModelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new LOTWWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        LOTWWallModelFvPatchScalarField
        (
            const LOTWWallModelFvPatchScalarField&,
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
                new LOTWWallModelFvPatchScalarField(*this, iF)
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
