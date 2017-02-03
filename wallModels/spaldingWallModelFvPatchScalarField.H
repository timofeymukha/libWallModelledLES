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
spaldingWallModel

Description
    Class for the wall model based on Spalding's law of the wall.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    spaldingWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef spaldingWallModelFvPatchScalarField_H
#define spaldingWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class spaldingWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class spaldingWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Data
        //RootFinder & rootFinder_;
    
        //- The "law of the wall" to be used
        autoPtr<LawOfTheWall> law_;
    
    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;
        
        //- Benchmark for calculating the friction velocity
        virtual tmp<scalarField> calcUTauBench(const scalarField& magGradU) const;


public:

    //- Runtime type information
    TypeName("spaldingWallModel");


    // Constructors

        //- Construct from patch and internal field
        spaldingWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaldingWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  spaldingWallModelFvPatchScalarField
        //  onto a new patch
        spaldingWallModelFvPatchScalarField
        (
            const spaldingWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaldingWallModelFvPatchScalarField
        (
            const spaldingWallModelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new spaldingWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        spaldingWallModelFvPatchScalarField
        (
            const spaldingWallModelFvPatchScalarField&,
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
                new spaldingWallModelFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif