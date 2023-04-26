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
    Foam::IntegratedReichardtLawOfTheWall

@brief
    Model based on integrating the law of the wall proposed by Reichardt.

    \f[
    u^+ = \frac{1}{\kappa } \ln (1 + \kappa y^+) + C \left( 1 - \exp (-y^+/B_1)-
          \frac{y^+}{B_1} \exp(-y^+/B_2) \right)
    \f]   

    Usage:
    \verbatim
    Law
    {
        type      IntegratedReichardt;
        kappa     value; (default 0.4)
        B1        value; (default 11)
        B2        value; (default 3)
        C         value; (default 7.8)
    }
    \endverbatim

    Reference:
    \verbatim
        Reichardt, H. (1951).
        Vollstandige Darstellung der turbulenten Geschwindigkeitsverteilung in
        glatten Leitungen.
        Zeitschrift fur Angewandte Mathematik und Mechanik 31(7) (pp. 208-219).
    \endverbatim

Contributors/Copyright:
    2018-2020 Timofey Mukha

SourceFiles
    IntegratedReichardtLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef IntegratedReichardtLawOfTheWall_H
#define IntegratedReichardtLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"
#include "MultiCellSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



/*---------------------------------------------------------------------------*\
                    Class IntegratedReichardtLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class IntegratedReichardtLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa_;
        
        //- The B1 model constant
        scalar B1_;
        
        //- The B1 model constant
        scalar B2_;
        
        //- The C model constant
        scalar C_;
        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("IntegratedReichardt");
#endif

    // Constructors

        //- Construct provided dictionary
        IntegratedReichardtLawOfTheWall
        (
            const dictionary &
        );
        
        //- Construct provided TypeName and dictionary
        IntegratedReichardtLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        IntegratedReichardtLawOfTheWall
        (
            const scalar kappa=0.4,
            const scalar B1=11,
            const scalar B2=3,
            const scalar C=7.8
        );

        // Destructor
        virtual ~IntegratedReichardtLawOfTheWall() {};
        
        //- Copy constructor
        IntegratedReichardtLawOfTheWall
        (
            const IntegratedReichardtLawOfTheWall &
        ) = default;

        //- Assignment
        IntegratedReichardtLawOfTheWall & operator=
        (
            const IntegratedReichardtLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const override
        {
            return autoPtr<LawOfTheWall>
            (
                new IntegratedReichardtLawOfTheWall(*this) 
            );
        }

    // Member Functions
        
        //- Return the kappa constant
        scalar kappa() const
        {
            return kappa_;
        }
        
        //- Return the B1 constant
        scalar B1() const
        {
            return B1_;
        }
        
        //- Return the B2 constant
        scalar B2() const
        {
            return B2_;
        }
        
        //- Return the C constant
        scalar C() const
        {
            return C_;
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

        scalar valueMulticell
        (
            const MultiCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const;

        scalar value
        (
            scalar uIntegral,
            scalar h1,
            scalar h2,
            scalar uTau,
            scalar nu
        ) const;

        //- Return the value of the  derivative of the implicit function
        //  defining the law.
        virtual scalar
        derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar derivativeMulticell
        (
            const MultiCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const;

        scalar
        derivative
        (
            scalar h1,
            scalar h2,
            scalar uTau,
            scalar nu
        ) const;

        //- The log-term in the integrated law
        scalar logTerm(scalar y, scalar uTau, scalar nu) const;
        
        //- The exp-term in the integrated law
        scalar expTerm(scalar y, scalar uTau, scalar nu) const;
        
        //- The derivative of the log-term in the integrated law wrt uTau
        scalar logTermDerivative(scalar y, scalar uTau, scalar nu) const;
        
        //- The derivative of the exp-term in the integrated law wrt uTau
        scalar expTermDerivative(scalar y, scalar uTau, scalar nu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
