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
    Foam::LSQRWallModelFvPatchScalarField

Description

Contributors/Copyright:
    2019 Timofey Mukha

SourceFiles
    LSQRWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef LSQRWallModelFvPatchScalarField_H
#define LSQRWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"
#include "LawOfTheWall.H"
#include "RootFinder.H"
#include "MultiCellSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class LSQRWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class LSQRWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:

    // Protected Data

        //- Pointer to the root finder
        autoPtr<RootFinder> rootFinder_;

        //- Pointer to the LSQR to be used
        autoPtr<LawOfTheWall> law_;

        //- Pointer to the sampler
        autoPtr<MultiCellSampler> sampler_;

    // Protected Member Functions
        //- Write root finder and LSQR properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        //- Create fields and add to registry
        virtual void createFields() const;
        
        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField & magGradU) const;
        

        
public:

    //- Runtime type information
    TypeName("LSQRWallModel");


    // Constructors

        //- Construct from patch and internal field
        LSQRWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        LSQRWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  LSQRWallModelFvPatchScalarField
        //  onto a new patch
        LSQRWallModelFvPatchScalarField
        (
            const LSQRWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        LSQRWallModelFvPatchScalarField
        (
            const LSQRWallModelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new LSQRWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        LSQRWallModelFvPatchScalarField
        (
            const LSQRWallModelFvPatchScalarField&,
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
                new LSQRWallModelFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        //- Write to stream
        virtual void write(Ostream& os) const;

        virtual void updateCoeffs();

        //- Return the sampler
        const MultiCellSampler & sampler() const
        {
            return sampler_();
        }

        MultiCellSampler & sampler()
        {
            return sampler_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
