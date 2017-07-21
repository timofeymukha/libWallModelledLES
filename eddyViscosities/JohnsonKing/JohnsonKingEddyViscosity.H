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
    Foam::JohnsonKingEddyViscosity

Description
    Eddy viscosity model based on the work of Johnson and King.
 
    Corresponds to the formulation used by Cabot in [Cabot].

    \f[
        \nu_t = \kappa u_\tau y (1 - \exp(- y^+/A))
    \f]

    Usage:
    \verbatim
    EddyViscosity
    {
        type      JohnsonKing;
        APlus     value; (default 19)
        kappa     value; (default 0.4)
    }
    \endverbatim

    References:
    \verbatim
        [JK]
        Johnson, D. A., & King, L. S. (1985).
        A mathematically simple turbulence closure model for attached and
        separated turbulent boundary layers.
        AIAA journal,
        23(11), 1684-1692.

        [Cabot]
        Cabot, W. (1995).
        Large-eddy simulations with wall models.
        Center for turbulence research,
        Annual Research Briefs, 1995
    \endverbatim

Authors
    Timofey Mukha, Saleh Rezaeiravesh.

SourceFiles
    JohnsonKingEddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonKingEddyViscosity_H
#define JohnsonKingEddyViscosity_H

#include "EddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class JohnsonKingEddyViscosity Declaration
\*---------------------------------------------------------------------------*/

class JohnsonKingEddyViscosity: public EddyViscosity
{
    
    // Private data
        
        //- The A+ model constant
        scalar APlus_;
        
        //- The kappa model constant
        scalar kappa_;
        
public:
    
    TypeName("JohnsonKing");

    // Constructors
        //- Construct with default values
    //    JohnsonKingEddyViscosity();
    
        //- Construct provided kappa, B and viscosity
      //  JohnsonKingEddyViscosity(scalar A, scalar B);
        
        //- Construct provided dictionary
        JohnsonKingEddyViscosity(const dictionary &);
        
        JohnsonKingEddyViscosity(const word &, const dictionary &);

    // Destructor
        virtual ~JohnsonKingEddyViscosity() {};

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
        value(const scalarList & y, scalar uTau, scalar nu) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
