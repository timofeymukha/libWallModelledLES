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
ODEWallModel

Description
    Class for wall models based on an ODE.

Authors
    Timofey Mukha.  All rights reserved.

SourceFiles
    ODEWallModelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ODEWallModelFvPatchScalarField_H
#define ODEWallModelFvPatchScalarField_H

#include "wallModelFvPatchScalarField.H"
#include "EddyViscosity.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class ODEWallModelPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ODEWallModelFvPatchScalarField
:
    public wallModelFvPatchScalarField
{
protected:
    
    // Protected Data
        //- Pointer to an eddy viscosity model    
        autoPtr<EddyViscosity> eddyViscosity_;
    
        //- A list of lists of points each defining a 1d mesh between on and h
        //  for a given face
        scalarListList meshes_;
        
        //- Maximum amount of iterations for coupling between uTau and nut
        label maxIter_;
        
        //- Error for exiting the uTau and nut coupling loop
        scalar eps_;
 
        //- number of points in the ODE mesh covering the distance h normal to the wall
        label nMeshY_;


    
    // Protected Member Functions
        //- Write root finder and ODE properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau() const;
        
        //- Create a 1d mesh for each face
        void createMeshes();
        
        //- Numerical integration
        scalar integrate(const scalarList & y, const scalarList & v) const;
        
        //source term
        virtual void source(const int, const vector, vector &) const
        {};

        
public:

    //- Runtime type information
    TypeName("ODEWallModel");


    // Constructors

        //- Construct from patch and internal field
        ODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ODEWallModelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ODEWallModelFvPatchScalarField
        //  onto a new patch
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );



    // Member functions

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
            
        scalar eps() const
        {
            return eps_;
        }

        scalar maxIter() const
        {
            return maxIter_;
        }

        scalar nMeshY() const
        {
            return nMeshY_;
        }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
