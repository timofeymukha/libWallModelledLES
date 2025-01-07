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
    Foam::SpaldingExplicitLawOfTheWall

@brief
    An explicit approximation of Spalding's law.

    \f[
    y^+ = u^+ +  \exp(-\kappa B) \left( \exp (\kappa B) - 1 - \kappa u^+ -
             \frac{1}{2}(\kappa u^+)^2 - \frac{1}{6}(\kappa u^+)^3 \right)
    \f]

    The approximation is constructed as a modification to the law of Cai
    and Sagaut.

    Usage:
    \verbatim
    Law
    {
        type      ExplicitSpalding;
        kappa     value; (default 0.4)
        B         value; (default 5)
    }
    \endverbatim

    Reference for the original law:
    \verbatim
        Spalding, D. B. (1961).
        A single formula for the 'law of the wall'.
        Journal of Applied Mechanics,
        28(3), 455-458.
    \endverbatim

Contributors/Copyright:
    2024 Timofey Mukha


SourceFiles
    SpaldingExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef SpaldingExplicitLawOfTheWall_H
#define SpaldingExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"
#include "CaiSagautExplicitLawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SpaldingExplicitLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class SpaldingExplicitLawOfTheWall: public ExplicitLawOfTheWall
{

    // Private data

        //- The kappa model constant
        scalar kappa_;

        //- The B model constant
        scalar B_;

        CaiSagautExplicitLawOfTheWall CaiSagaut_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Spalding");
#endif

    // Constructors

        //- Construct provided dictionary
        SpaldingExplicitLawOfTheWall(const dictionary &);

        //- Construct provided TypeName and dictionary
        SpaldingExplicitLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        SpaldingExplicitLawOfTheWall
        (
            const scalar kappa,
            const scalar B
        );


        //- Copy constructor
        SpaldingExplicitLawOfTheWall
        (
            const SpaldingExplicitLawOfTheWall &
        ) = default;

        //- Assignment
        SpaldingExplicitLawOfTheWall & operator=
        (
            const SpaldingExplicitLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<ExplicitLawOfTheWall> clone() const override
        {
            return autoPtr<ExplicitLawOfTheWall>
            (
                new SpaldingExplicitLawOfTheWall(*this)
            );
        }

    // Destructor
        virtual ~SpaldingExplicitLawOfTheWall() {};

    // Member Functions

        //- Return the kappa constant
        scalar kappa() const
        {
            return kappa_;
        }

        //- Return the B constant
        scalar B() const
        {
            return B_;
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
