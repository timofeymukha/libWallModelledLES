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
    Foam::SpaldingLawOfTheWall

@brief
    The law of the wall proposed by Spalding.

    \f[
    y^+ = u^+ +  \exp(-\kappa B) \left( \exp (\kappa B) - 1 - \kappa u^+ -
             \frac{1}{2}(\kappa u^+)^2 - \frac{1}{6}(\kappa u^+)^3 \right)
    \f]

    Since the formula above is differentiable, it is advised to use Newton's
    method to solve for \f$u_\tau \f$.

    Usage:
    \verbatim
    Law
    {
        type      Spalding;
        kappa     value; (default 0.4)
        B         value; (default 5.5)
    }
    \endverbatim
    
    Reference:
    \verbatim
        Spalding, D. B. (1961).
        A single formula for the 'law of the wall'.
        Journal of Applied Mechanics,
        28(3), 455-458.
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh


SourceFiles
    SpaldingLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef SpaldingLawOfTheWall_H
#define SpaldingLawOfTheWall_H

#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SpaldingLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class SpaldingLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa_;
        
        //- The B model constant
        scalar B_;
        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Spalding");
#endif

    // Constructors
      
        //- Construct provided dictionary
        SpaldingLawOfTheWall(const dictionary &);
        
        //- Construct provided TypeName and dictionary
        SpaldingLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        SpaldingLawOfTheWall
        (
            const scalar kappa,
            const scalar B
        );

        //- Copy constructor
        SpaldingLawOfTheWall(const SpaldingLawOfTheWall &) = default;

        //- Assignment
        SpaldingLawOfTheWall & operator=(const SpaldingLawOfTheWall &) = default;

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const override
        {
            return autoPtr<LawOfTheWall>
            (
                new SpaldingLawOfTheWall(*this) 
            );
        }

    //- Destructor
        virtual ~SpaldingLawOfTheWall() {};

    // Member Functions
        
        scalar kappa() const
        {
            return kappa_;
        }
        
        scalar B() const
        {
            return B_;
        }
        
        //- Print the model coefficients to terminal
        virtual void printCoeffs() const override;

        //- Return the value of the implicit function defining the law
        virtual scalar value
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar value(scalar u, scalar y, scalar uTau, scalar nu) const;

                             
        //- Return the value of the  derivative of the implicit function
        //  defining the law.
        virtual scalar derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar derivative(scalar u, scalar y, scalar uTau, scalar nu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
