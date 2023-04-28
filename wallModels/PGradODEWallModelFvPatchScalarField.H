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
    Foam::PGradODEWallModelFvPatchScalarField

@brief
    ODE wall model with the source term equal to the wall-parallel pressure
    gradient.

    Usage
    \verbatim
    patchName
    {
        type                PGradODEWallModel;
        value               uniform 0;
        maxIter             value; (default 10)
        eps                 value; (default 1e-3)
        nMeshY              value; (default 30)

        EddyViscosity 
        {
            type            EddyViscosityType;
            otherParams     value;
        }
    }
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    PGradODEWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef PGradODEWallModelFvPatchScalarField_H
#define PGradODEWallModelFvPatchScalarField_H

#include "ODEWallModelFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class PGradODEWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class PGradODEWallModelFvPatchScalarField
:
    public ODEWallModelFvPatchScalarField
{
protected:

    // Protected Data
        
    // Protected Member Functions
        
        //- Write root finder and PGradODE properties to stream
        virtual void writeLocalEntries(Ostream &) const;
              
        //- Compute the source term based on the pressure gradient
        virtual void source(vectorField &) const;
        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("PGradODEWallModel");
#endif


    // Constructors

        //- Construct from patch and internal field
        PGradODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        PGradODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  PGradODEWallModelFvPatchScalarField
        //  onto a new patch
        PGradODEWallModelFvPatchScalarField
        (
            const PGradODEWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        PGradODEWallModelFvPatchScalarField
        (
            const PGradODEWallModelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new PGradODEWallModelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        PGradODEWallModelFvPatchScalarField
        (
            const PGradODEWallModelFvPatchScalarField&,
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
                new PGradODEWallModelFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

    virtual void write(Ostream& os) const;   
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
