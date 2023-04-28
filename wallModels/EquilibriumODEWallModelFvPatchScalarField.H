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
    Foam::EquilibriumODEWallModelFvPatchScalarField

@brief
    ODE model with a zero source-term.
  
    Physically, this means that the total shear stress is assumed to be constant
    from 0 to h. In practice, this model is equivalent to using an equilibrium
    algebraic model, but using numerical integration intead of Newton's solver
    to get the shear stress. The former is numerically stable, so this model
    probably never crash. This is convenient, but also makes monitoring your
    results as the simulation goes more important.

    Usage
    \verbatim
    patchName
    {
        type                EquilibriumODEWallModel;
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
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh

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
        //- Write model properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- The source term, which is zero
        virtual void source(vectorField &) const
        {};
        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("EquilibriumODEWallModel");
#endif


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
        
        //- Write
        virtual void write(Ostream& os) const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
