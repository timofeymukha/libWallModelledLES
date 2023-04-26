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
    Foam::LOTWWallModelFvPatchScalarField

@brief
    Class for the wall models based on a law of the wall.
 
    Needs a root finder and a law of the wall to be defined.

    Usage
    \verbatim
    patchName
    {
        type                LOTWWallModel;
        value               uniform 0;
        RootFinder
        {
            type            RootFinderType;
            otherParams     value;
        }

        Law
        {
            type            LawOfTheWallType;
            otherParams     value;
        }
    }
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    LOTWWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef LOTWWallModelFvPatchScalarField_H
#define LOTWWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class RootFinder;
class LawOfTheWall;
class SingleCellSampler;

/*---------------------------------------------------------------------------*\
          Class LOTWWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class LOTWWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Data
    
        //- Pointer to the root finder
        autoPtr<RootFinder> rootFinder_;
    
        //- Pointer to the LOTW to be used
        autoPtr<LawOfTheWall> law_;

        //- The sampler
        autoPtr<SingleCellSampler> sampler_;
    
    // Protected Member Functions
        //- Write root finder and LOTW properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField & magGradU) const;
        

        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("LOTWWallModel");
#endif


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
                
        #ifdef FOAM_FVPATCHFIELD_NO_COPY
        #else
        //- Construct as copy
        LOTWWallModelFvPatchScalarField
        (
            const LOTWWallModelFvPatchScalarField&
        );
        #endif

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
