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
    Foam::KnownWallShearStressWallModelFvPatchScalarField

@brief
    Wall model that uses the values of the magnitude of the mean wall shear
    stress known a priori, for example, from a RANS precursor or an analytical
    solution.
    
    Requires field tauWall with the values of the magnitude of the mean wall
    shear stress to be present.

Contributors/Copyright:
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    KnownWallShearStressWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef KnownWallShearStressWallModelFvPatchScalarField_H
#define KnownWallShearStressWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class KnownWallShearStressWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class KnownWallShearStressWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Member Functions
        //- Write model properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity at the wall
        virtual tmp<scalarField> calcNut() const;
               
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("KnownWallShearStressWallModel");
#endif

    // Constructors

        //- Construct from patch and internal field
        KnownWallShearStressWallModelFvPatchScalarField
        (
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &
        );

        //- Construct from patch, internal field and dictionary
        KnownWallShearStressWallModelFvPatchScalarField
        (
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &
        );

        //- Construct by mapping given KnownWallShearStressWallModelFvPatchScalarField
        //  onto a new patch
        KnownWallShearStressWallModelFvPatchScalarField
        (
            const KnownWallShearStressWallModelFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &
        );

        //- Construct as copy
        KnownWallShearStressWallModelFvPatchScalarField
        (
            const KnownWallShearStressWallModelFvPatchScalarField &
        );

        //- Construct as copy setting internal field reference
        KnownWallShearStressWallModelFvPatchScalarField
        (
            const KnownWallShearStressWallModelFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &
        );
        

    // Member functions


        //- Write
        virtual void write(Ostream& os) const;            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
