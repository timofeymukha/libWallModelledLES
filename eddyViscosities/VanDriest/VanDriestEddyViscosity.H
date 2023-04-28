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

@brief
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

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

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
        
        //- The kappa model constant
        scalar kappa_;

        //- The A+ model constant
        scalar APlus_;

        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("VanDriest");
#endif

    // Constructors
      
        //- Construct provided dictionary and sampler
        VanDriestEddyViscosity(const dictionary &);
        
        //- Construct provided name, dictionary and sampler
        VanDriestEddyViscosity(const word &, const dictionary &);

        //- Construct from model cosntants
        VanDriestEddyViscosity(const scalar kappa, const scalar APlus);

        //- Copy constructor
        VanDriestEddyViscosity(const VanDriestEddyViscosity &) = default;

        //- Assignment
        VanDriestEddyViscosity & operator=
        (
            const VanDriestEddyViscosity &
        ) = default;

        //- Clone
        virtual autoPtr<EddyViscosity> clone() const override
        {
            return autoPtr<EddyViscosity>
            (
                new VanDriestEddyViscosity(*this)
            );
        }

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
        virtual void printCoeffs() const override;

        //- Return the value of nut
        virtual scalarList
        value
        (
            const SingleCellSampler & sampler,
            const label index,
            const scalarList & y,
            const scalar uTau,
            const scalar nu
        ) const override;

        scalarList value
        (
            const scalarList & y,
            const scalar uTau,
            const scalar nu
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
