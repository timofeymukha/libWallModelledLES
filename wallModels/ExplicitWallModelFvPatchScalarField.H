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
    Foam::ExplicitWallModelFvPatchScalarField

@brief
    Class for the wall models based on a an explicit expression for the wall
    stress.

    Needs an explicit expression of the wall to be defined, which we refer to
    as the "law" in similarity with algebraic models.

    Usage
    \verbatim
    patchName
    {
        type                ExpliciWallModel;
        value               uniform 0;

        Law
        {
            type            ExpressionType;
            otherParams     value;
        }
    }
    \endverbatim

Contributors/Copyright:
    2024-2025 Timofey Mukha

SourceFiles
    ExplicitWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitWallModelFvPatchScalarField_H
#define ExplicitWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"
#include <boost/math/tools/roots.hpp>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class RootFinder;
class ExplicitLawOfTheWall;
class SingleCellSampler;


/*---------------------------------------------------------------------------*\
          Class ExplicitWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ExplicitWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:

    // Protected Data

        //- Pointer to the ExplicitLOTW to be used
        autoPtr<ExplicitLawOfTheWall> law_;

        //- The sampler
        autoPtr<SingleCellSampler> sampler_;

    // Protected Member Functions
        //- Write root finder and ExplicitLOTW properties to stream
        virtual void writeLocalEntries(Ostream &) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField & magGradU) const;



public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("ExplicitWallModel");
#endif


    // Constructors

        //- Construct from patch and internal field
        ExplicitWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ExplicitWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ExplicitWallModelFvPatchScalarField
        //  onto a new patch
        ExplicitWallModelFvPatchScalarField
        (
            const ExplicitWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        #ifdef FOAM_FVPATCHFIELD_NO_COPY
        #else
        //- Construct as copy
        ExplicitWallModelFvPatchScalarField
        (
            const ExplicitWallModelFvPatchScalarField&
        );
        #endif

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ExplicitWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ExplicitWallModelFvPatchScalarField
        (
            const ExplicitWallModelFvPatchScalarField&,
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
                new ExplicitWallModelFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        SingleCellSampler & sampler()
        {
            return sampler_();
        }

        const SingleCellSampler & sampler() const
        {
            return sampler_();
        }

        virtual void updateCoeffs();

        //- Write to stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
