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
    Foam::VanDriestEddyViscosity

Description
    A mixing-length based eddy viscosity model using a van Driest damping
    function.
    Corresponds to the formulation used by Cabot in [Cabot].

    \f[
        \nu_t = \kappa u_\tau y (1 - \exp(- y^+/A))
    \f]

    Usage:
    \verbatim
    EddyViscosity
    {
        type      VanDriest;
        APlus     value; (default 18)
        kappa     value; (default 0.4)
    }
    \endverbatim

    References:
    \verbatim
        [Cabot]
        Cabot, W. (1995).
        Large-eddy simulations with wall models.
        Center for turbulence research,
        Annual Research Briefs, 1995
    \endverbatim

Authors
    Timofey Mukha, Saleh Rezaeiravesh.

SourceFiles
    VanDriestEddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef VanDriestEddyViscosity_H
#define VanDriestEddyViscosity_H

#include "EddyViscosity.H"
#include "Sampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class VanDriestEddyViscosity Declaration
\*---------------------------------------------------------------------------*/

class VanDriestEddyViscosity: public EddyViscosity
{
    
    // Private data
        
        //- The A+ model constant
        scalar APlus_;
        
        //- The kappa model constant
        scalar kappa_;
        
public:
    
    TypeName("VanDriest");

    // Constructors
      
        //- Construct provided dictionary and sampler
        VanDriestEddyViscosity(const dictionary &, Sampler &);
        
        //- Construct provided name, dictionary and sampler
        VanDriestEddyViscosity(const word &, const dictionary &, Sampler &);

    // Destructor
        virtual ~VanDriestEddyViscosity() {};

    // Member Functions
        
        scalar APlus() const
        {
            return APlus_;
        }
        
        scalar kappa() const
        {
            return kappa_;
        }
        
        //- Print info
        virtual void printCoeffs() const;

        //- Return the value of nut
        virtual scalarList
        value(label index, const scalarList & y, scalar uTau, scalar nu) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
