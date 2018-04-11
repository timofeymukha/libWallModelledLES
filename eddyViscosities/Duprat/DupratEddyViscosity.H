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

Description
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

Authors
    Timofey Mukha, Saleh Rezaeiravesh.

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
        
        //- The A+ model constant
        scalar APlus_;
        
        //- The kappa model constant
        scalar kappa_;
        
        //- The beta model cosntant
        scalar beta_;
        
public:
    
    TypeName("Duprat");

    // Constructors
        
        //- Construct provided dictionary and sampler
        DupratEddyViscosity(const dictionary &, Sampler &);
        
        //- Construct provided name, dictionary and sampler
        DupratEddyViscosity(const word &, const dictionary &, Sampler &);

    // Destructor
        virtual ~DupratEddyViscosity() {};

    // Member Functions
        
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
        virtual void printCoeffs() const;

        //- Return the value of nut
        virtual scalarList
        value(label index, const scalarList & y, scalar uTau, scalar nu) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif