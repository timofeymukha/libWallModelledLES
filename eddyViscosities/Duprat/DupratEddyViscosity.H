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
    Foam::DupratEddyViscosity

@brief
    Eddy viscosity model based on the work of Duprat et al. [Duprat].

    \f[
        \nu_t = \nu \kappa y^* [\alpha + y^*(1 - \alpha)^{3/2}]^\beta 
                (1 - \exp(-y^* /(1 + A \alpha^3)))^2
    \f]

    Usage:
    \verbatim
    EddyViscosity
    {
        type      Duprat;
        APlus     value; (default 17)
        beta      value; (default 0.78)
        kappa     value; (default 0.4)
    }
    \endverbatim

    References:
    \verbatim
        [Duprat]
        Duprat, C., Balarac, G., Meetais, O.,
        Congedo, P. M., & Brugiere, O. (2011).
        A wall-layer model for large-eddy simulations of turbulent flows
        with/out pressure gradient.
        Physics of Fluids, 
        23(1), 15101.
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    DupratEddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef DupratEddyViscosity_H
#define DupratEddyViscosity_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "EddyViscosity.H"
#include "Sampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DupratEddyViscosity Declaration
\*---------------------------------------------------------------------------*/

class DupratEddyViscosity: public EddyViscosity
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa_;

        //- The A+ model constant
        scalar APlus_;
               
        //- The beta model cosntant
        scalar beta_;
        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Duprat");
#endif

    // Constructors
        
        //- Construct provided dictionary and sampler
        DupratEddyViscosity(const dictionary &);
        
        //- Construct provided name, dictionary and sampler
        DupratEddyViscosity(const word &, const dictionary &);

        //- Construct from model constants
        DupratEddyViscosity
        (
            const scalar kappa,
            const scalar APlus,
            const scalar beta
        );

        //- Copy constructor
        DupratEddyViscosity(const DupratEddyViscosity &) = default;

        //- Assignment
        DupratEddyViscosity & operator=(const DupratEddyViscosity &) = default;

        //- Clone
        virtual autoPtr<EddyViscosity> clone() const override
        {
            return autoPtr<EddyViscosity>
            (
                new DupratEddyViscosity(*this)
            );
        }

    // Destructor
        virtual ~DupratEddyViscosity() {};

    // Member Functions

        //- Add necessary sampled fields to the sampler
        virtual void addFieldsToSampler(Sampler & sampler) const override;
        
        scalar APlus() const
        {
            return APlus_;
        }
        
        scalar kappa() const
        {
            return kappa_;
        }
        
        scalar beta() const
        {
            return beta_;
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
            const scalar magPGrad,
            const scalar uTau,
            const scalar nu
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
