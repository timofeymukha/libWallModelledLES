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
    Foam::WernerWengleLawOfTheWall

Description
    The law of the wall proposed by Werner and Wengle.

    \f{eqnarray*}{
        u^+ &=& y^+ \quad y+ \leq 11.81\\
        u^+ &=& A(y^+)^B \quad y+ > 11.81

    \f}    

    Usage:
    \verbatim
    Law
    {
        type      WernerWengle;
        A         value; (default 8.3)
        B         value; (default 1/7)
    }
    \endverbatim

    Reference:
    \verbatim
        Werner, H., & Wengle, H. (1993).
        Large-eddy simulation of turbulent flow over and around a cube in a
        plate channel.
        In Turbulent Shear Flows 8 (pp. 155-168).
        Springer, Berlin, Heidelberg.
    \endverbatim

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    WernerWengleLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef WernerWengleLawOfTheWall_H
#define WernerWengleLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class WernerWengleLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class WernerWengleLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar A_;
        
        //- The B model constant
        scalar B_;
        
public:
    
    TypeName("WernerWengle");

    // Constructors

        //- Construct provided dictionary and cell index list
        WernerWengleLawOfTheWall(const dictionary &, const Foam::CellIndexList &);
        
        //- Construct provided TypeName and dictionary and cell index list
        WernerWengleLawOfTheWall
        ( 
            const word & lawname,
            const dictionary &,
            const Foam::CellIndexList &
        );

    // Destructor
        virtual ~WernerWengleLawOfTheWall() {};

    // Member Functions
        
        //- Return the A constant
        scalar A() const
        {
            return A_;
        }
        
        //- Return the B constant
        scalar B() const
        {
            return B_;
        }
        
        //- Print info to terminal
        virtual void printCoeffs() const;

        //- Return the value of the implicit function defining the law
        virtual scalar
        value(scalar u, scalar index, scalar uTau, scalar nu) const;
        
        //- Return the value of the derivative of the implicit function
        //  defining the law
        virtual scalar
        derivative(scalar u, scalar index, scalar uTau, scalar nu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
