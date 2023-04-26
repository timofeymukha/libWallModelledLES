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
    Foam::IntegratedWernerWengleLawOfTheWall

@brief
    The integrated version of the wall proposed by Werner and Wengle.
    Explicitely couples shear stress with the wall-normal-averaged velocity.
    Derived assuming that the sampling is from the wall-adjacent cell.

    Based on the following law of the wall.

    \f{eqnarray*}{
        u^+ &=& y^+ \quad y+ \leq 11.81\\
        u^+ &=& A(y^+)^B \quad y^+ > 11.81
    \f}    

    Usage:
    \verbatim
    Law
    {
        type      IntegratedIntegratedWernerWengle;
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

Contributors/Copyright:
    2018-2019 Timofey Mukha

SourceFiles
    IntegratedWernerWengleLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef IntegratedWernerWengleLawOfTheWall_H
#define IntegratedWernerWengleLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class IntegratedWernerWengleLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class IntegratedWernerWengleLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar A_;
        
        //- The B model constant
        scalar B_;
        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("IntegratedWernerWengle");
#endif

    // Constructors

        //- Construct provided dictionary
        IntegratedWernerWengleLawOfTheWall
        (
            const dictionary &
        );
        
        //- Construct provided TypeName and dictionary
        IntegratedWernerWengleLawOfTheWall
        ( 
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        IntegratedWernerWengleLawOfTheWall
        ( 
            const scalar A,
            const scalar B
        );

        //- Copy constructor
        IntegratedWernerWengleLawOfTheWall
        (
            const IntegratedWernerWengleLawOfTheWall &
        ) = default;

        //- Assignment
        IntegratedWernerWengleLawOfTheWall & operator=
        (
            const IntegratedWernerWengleLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const override
        {
            return autoPtr<LawOfTheWall>
            (
                new IntegratedWernerWengleLawOfTheWall(*this) 
            );
        }

        //- Destructor
        virtual ~IntegratedWernerWengleLawOfTheWall() {};

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
        virtual void printCoeffs() const override;

        //- Return the value of the implicit function defining the law
        virtual scalar value
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar value
        (
            scalar u,
            scalar y,
            scalar uTau,
            scalar nu
        ) const;

        //- Return the value of the  derivative of the implicit function
        //  defining the law.
        virtual scalar derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar derivative() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
