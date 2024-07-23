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
    Foam::ReichardtExplicitLawOfTheWall

@brief
    An explicit approximation of Reichardt's law.

    \f[
    u^+ = \frac{1}{\kappa } \ln (1 + \kappa y^+) + C \left( 1 - \exp (-y^+/B_1)-
          \frac{y^+}{B_1} \exp(-y^+/B_2) \right)
    \f]

    The approximation is constructed as a modification to the law of Cai
    and Sagaut.

    Usage:
    \verbatim
    Law
    {
        type      ExplicitReichardt;
        kappa     value; (default 0.4)
        B1        value; (default 11)
        B2        value; (default 3)
        C         value; (default 7.8)
    }
    \endverbatim

    Reference for the original law:
    \verbatim
        Reichardt, H. (1951).
        Vollstandige Darstellung der turbulenten Geschwindigkeitsverteilung in
        glatten Leitungen.
        Zeitschrift fur Angewandte Mathematik und Mechanik 31(7) (pp. 208-219).
    \endverbatim

Contributors/Copyright:
    2024 Timofey Mukha


SourceFiles
    ReichardtExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef ReichardtExplicitLawOfTheWall_H
#define ReichardtExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"
#include "CaiSagautExplicitLawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ReichardtExplicitLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class ReichardtExplicitLawOfTheWall: public ExplicitLawOfTheWall
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

        CaiSagautExplicitLawOfTheWall CaiSagaut_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Reichardt");
#endif

    // Constructors

        //- Construct provided dictionary
        ReichardtExplicitLawOfTheWall(const dictionary &);

        //- Construct provided TypeName and dictionary
        ReichardtExplicitLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        ReichardtExplicitLawOfTheWall
        (
            const scalar kappa,
            const scalar B1,
            const scalar B2,
            const scalar C
        );


        //- Copy constructor
        ReichardtExplicitLawOfTheWall
        (
            const ReichardtExplicitLawOfTheWall &
        ) = default;

        //- Assignment
        ReichardtExplicitLawOfTheWall & operator=
        (
            const ReichardtExplicitLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<ExplicitLawOfTheWall> clone() const override
        {
            return autoPtr<ExplicitLawOfTheWall>
            (
                new ReichardtExplicitLawOfTheWall(*this)
            );
        }

    // Destructor
        virtual ~ReichardtExplicitLawOfTheWall() {};

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
        virtual scalar uTau
        (
            const SingleCellSampler & sampler,
            label index,
            scalar nu
        ) const override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
