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
    Foam::CaiSagautExplicitLawOfTheWall

@brief
    The explicit law of the wall proposed by Cai and Sagaut. This is "model 2"
    in the paper. Uses the Lmabert W function to explicitly represent the
    log-law.

    \f[
    u^+ = \exp(-Re_y / s)^p \sqrt(Re_y) +
          (1 - \exp(-Re_y /s))^p W(Re_y \kappa E) / \kappa
    \f]

    Usage:
    \verbatim
    Law
    {
        type      CaiSagaut;
        kappa     value; (default 0.4)
        B         value; (default 5.5)
        p         value; (default 1.138)
        s         value; (default 217.8)
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
    CaiSagautExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef CaiSagautExplicitLawOfTheWall_H
#define CaiSagautExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CaiSagautExplicitLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class CaiSagautExplicitLawOfTheWall: public ExplicitLawOfTheWall
{

    // Private data

        //- The kappa model constant
        scalar kappa_;

        //- The B model constant
        scalar B_;

        //- The p model constant
        scalar p_;

        //- The s model constant
        scalar s_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("CaiSagaut");
#endif

    // Constructors

        //- Construct provided dictionary
        CaiSagautExplicitLawOfTheWall(const dictionary &);

        //- Construct provided TypeName and dictionary
        CaiSagautExplicitLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        CaiSagautExplicitLawOfTheWall
        (
            const scalar kappa,
            const scalar B,
            const scalar p,
            const scalar s
        );


        //- Copy constructor
        CaiSagautExplicitLawOfTheWall
        (
            const CaiSagautExplicitLawOfTheWall &
        ) = default;

        //- Assignment
        CaiSagautExplicitLawOfTheWall & operator=
        (
            const CaiSagautExplicitLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<ExplicitLawOfTheWall> clone() const override
        {
            return autoPtr<ExplicitLawOfTheWall>
            (
                new CaiSagautExplicitLawOfTheWall(*this)
            );
        }

    // Destructor
        virtual ~CaiSagautExplicitLawOfTheWall() {};

    // Member Functions

        // Return kappa
        scalar kappa() const
        {
            return kappa_;
        }

        // Return B
        scalar B() const
        {
            return B_;
        }

        // Return j
        scalar p() const
        {
            return p_;
        }

        // Return s
        scalar s() const
        {
            return s_;
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
