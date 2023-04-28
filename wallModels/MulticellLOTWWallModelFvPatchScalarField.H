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
    Foam::MulticellLOTWWallModelFvPatchScalarField

@brief
    Class for the wall models based on a law of the wall, sampling data
    from multiple cells. Currently hardcoded to use the integrated
    Reichardt law.
 
    Needs a root finder to be defined.

    Usage
    \verbatim
    patchName
    {
        type                MulticellLOTWWallModel;
        value               uniform 0;
        RootFinder
        {
            type            RootFinderType;
            otherParams     value;
        }
    }
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    MulticellLOTWWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef MulticellLOTWWallModelFvPatchScalarField_H
#define MulticellLOTWWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class RootFinder;
class IntegratedReichardtLawOfTheWall;
class MultiCellSampler;

/*---------------------------------------------------------------------------*\
          Class MulticellLOTWWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class MulticellLOTWWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Data
    
        //- Pointer to the root finder
        autoPtr<RootFinder> rootFinder_;
    
        //- Pointer to the MulticellLOTW to be used
        autoPtr<IntegratedReichardtLawOfTheWall> law_;

        //- The sampler
        autoPtr<MultiCellSampler> sampler_;
    
    // Protected Member Functions
        //- Write root finder and MulticellLOTW properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField & magGradU) const;
        

        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("MulticellLOTWWallModel");
#endif


    // Constructors

        //- Construct from patch and internal field
        MulticellLOTWWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        MulticellLOTWWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  MulticellLOTWWallModelFvPatchScalarField
        //  onto a new patch
        MulticellLOTWWallModelFvPatchScalarField
        (
            const MulticellLOTWWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
                
        #ifdef FOAM_FVPATCHFIELD_NO_COPY
        #else
        //- Construct as copy
        MulticellLOTWWallModelFvPatchScalarField
        (
            const MulticellLOTWWallModelFvPatchScalarField&
        );
        #endif

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new MulticellLOTWWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        MulticellLOTWWallModelFvPatchScalarField
        (
            const MulticellLOTWWallModelFvPatchScalarField&,
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
                new MulticellLOTWWallModelFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        MultiCellSampler & sampler()
        {
            return sampler_();
        }

        const MultiCellSampler & sampler() const
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
