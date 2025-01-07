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
    Foam::EquilibriumODEExplicitLawOfTheWall

@brief
    An explicit approximation of Equilibruim ODE model.

Contributors/Copyright:
    2024 Timofey Mukha


SourceFiles
    EquilibriumODEExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef EquilibriumODEExplicitLawOfTheWall_H
#define EquilibriumODEExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"
#include "CaiSagautExplicitLawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class EquilibriumODEExplicitLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class EquilibriumODEExplicitLawOfTheWall: public ExplicitLawOfTheWall
{

    // Private data

        //- The kappa model constant
        scalar kappa_;

        //- The Aplus model constant
        scalar Aplus_;

        CaiSagautExplicitLawOfTheWall CaiSagaut_;

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("EquilibriumODE");
#endif

    // Constructors

        //- Construct provided dictionary
        EquilibriumODEExplicitLawOfTheWall(const dictionary &);

        //- Construct provided TypeName and dictionary
        EquilibriumODEExplicitLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        EquilibriumODEExplicitLawOfTheWall
        (
            const scalar kappa,
            const scalar Aplus
        );


        //- Copy constructor
        EquilibriumODEExplicitLawOfTheWall
        (
            const EquilibriumODEExplicitLawOfTheWall &
        ) = default;

        //- Assignment
        EquilibriumODEExplicitLawOfTheWall & operator=
        (
            const EquilibriumODEExplicitLawOfTheWall &
        ) = default;

        //- Clone
        virtual autoPtr<ExplicitLawOfTheWall> clone() const override
        {
            return autoPtr<ExplicitLawOfTheWall>
            (
                new EquilibriumODEExplicitLawOfTheWall(*this)
            );
        }

    // Destructor
        virtual ~EquilibriumODEExplicitLawOfTheWall() {};

    // Member Functions

        //- Return the kappa constant
        scalar kappa() const
        {
            return kappa_;
        }

        //- Return the Aplus constant
        scalar Aplus() const
        {
            return Aplus_;
        }


        //- Print info to terminal
        virtual void printCoeffs() const override;

        //- Return the value of the friction velocity
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
