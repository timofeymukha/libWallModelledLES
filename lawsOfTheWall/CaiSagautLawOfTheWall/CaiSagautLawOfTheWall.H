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
    Foam::CaiSagautLawOfTheWall

@brief
    The explicit law of the wall proposed by Cai and Sagaut. This is "model 2"
    in the paper.

    \f[
    u^+ = \frac{1}{\kappa } \ln (1 + \kappa y^+) + C \left( 1 - \exp (-y^+/B_1)-
          \frac{y^+}{B_1} \exp(-y^+/B_2) \right)
    \f]

    Usage:
    \verbatim
    Law
    {
        type      CaiSagaut;
        kappa     value; (default 0.4)
        B        value; (default 5.5)
    }
    \endverbatim

    Reference:
    \verbatim
        Cai and Sagaut. (2021).
        Explicit wall models for large eddy simulation
        Physics of Fluids, 33 (4), 41703
    \endverbatim

Contributors/Copyright:
    2024 Timofey Mukha


SourceFiles
    CaiSagautLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef CaiSagautLawOfTheWall_H
#define CaiSagautLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CaiSagautLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class CaiSagautLawOfTheWall: public LawOfTheWall
{

    // Private data

        //- The kappa model constant
        scalar kappa_;

        //- The B model constant
        scalar B_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("CaiSagaut");
#endif

    // Constructors

        //- Construct provided dictionary
        CaiSagautLawOfTheWall(const dictionary &);

        //- Construct provided TypeName and dictionary
        CaiSagautLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        CaiSagautLawOfTheWall
        (
            const scalar kappa,
            const scalar B
        );


        //- Copy constructor
        CaiSagautLawOfTheWall
        (
            const CaiSagautLawOfTheWall &
        ) = default;

        //- Assignment
        CaiSagautLawOfTheWall & operator=
        (
            const CaiSagautLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const override
        {
            return autoPtr<LawOfTheWall>
            (
                new CaiSagautLawOfTheWall(*this)
            );
        }

    // Destructor
        virtual ~CaiSagautLawOfTheWall() {};

    // Member Functions

        //- Return the kappa constant
        scalar kappa() const
        {
            return kappa_;
        }

        //- Return the B1 constant
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

        //- Return the value of the derivative of the implicit function
        //  defining the law
        virtual scalar derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar derivative
        (
            scalar u,
            scalar y,
            scalar uTau,
            scalar nu
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
