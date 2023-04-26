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
    Foam::ODEWallModelFvPatchScalarField

@brief
    Base abstract class for ODE-based wall models.
  
    All the models are based on the following equation

    \f[
        \frac{\partial}{\partial y}
        \left[ (\nu + \nu_t)\frac{\partial U_i}{\partial y}\right] = F_i,
    \f]

    where i corresponds to the two wall-parallel coordinates.
    If the source term F is not dependent on y this can be integrated to give
    \f[
         \tau_{w,i} = 
            \left(
                U_i|_h -
                F_i \int^h_0 \frac{y}{\nu + \nu_t}dy
            \right)
            \bigg  / 
            \int^h_0 \frac{dy}{\nu + \nu_t}
    \f]

    The magnitude of the wall shear stress can then be expressed through the
    magnitude of the the wall-parallel component of velocity and the source
    term.
    Note that no ODE is actually being solved, only numerical integration.

    All ODE wall models require an eddy viscosity model. All share the following
    parameters
    - maxIter, the amount of iterations in the coupling loop between the wall
    shear stress and the eddy viscosity values.
    - eps, the relative error tolerance for the convergence of the wall shear
    stress.
    - nMeshY, the amount of nodes in the 1D grid discretising the distance
    between 0 and h.

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

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

        //- The sampler
        autoPtr<SingleCellSampler> sampler_;

        //- A list of lists of points each defining a 1d mesh between 0 and h
        //  for a given face
        scalarListList meshes_;
        
        //- Maximum amount of iterations for coupling between uTau and nut
        label maxIter_;
        
        //- Error for exiting the uTau and nut coupling loop
        scalar eps_;
 
        //- number of points in each mesh in meshes_
        label nMeshY_;

   
    // Protected Member Functions
        //- Write model properties to stream
        virtual void writeLocalEntries(Ostream &) const;
        
        //- Calculate the turbulence viscosity at the wall
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField & magGradU) const;
        
        //- Create a 1d mesh for each face
        void createMeshes();
        
        //- Numerical integration
        scalar integrate(const scalarList & y, const scalarList & v) const;
        
        //- Source term defining the type of ODE model
        virtual void source(vectorField &) const = 0;
        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
    TypeName("ODEWallModel");
#endif

    // Constructors

        //- Construct from patch and internal field
        ODEWallModelFvPatchScalarField
        (
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &
        );

        //- Construct from patch, internal field and dictionary
        ODEWallModelFvPatchScalarField
        (
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &
        );

        //- Construct by mapping given ODEWallModelFvPatchScalarField
        //  onto a new patch
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &
        );
        
        #ifdef FOAM_FVPATCHFIELD_NO_COPY
        #else
        //- Construct as copy
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField &
        );
        #endif

        //- Construct as copy setting internal field reference
        ODEWallModelFvPatchScalarField
        (
            const ODEWallModelFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &
        );
        

    // Member functions


        //- Write
        virtual void write(Ostream& os) const;
            
        virtual void updateCoeffs();

        //- Return the error tolerance
        scalar eps() const
        {
            return eps_;
        }

        //- Return the max number of coupling iterations
        scalar maxIter() const
        {
            return maxIter_;
        }

        //- Return the size of the 1d meshes
        scalar nMeshY() const
        {
            return nMeshY_;
        }

        SingleCellSampler & sampler()
        {
            return sampler_();
        }

        const SingleCellSampler & sampler() const
        {
            return sampler_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
