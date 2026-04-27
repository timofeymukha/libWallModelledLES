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
    An explicit approximation of the equilibrium ODE wall model.

    The equilibrium ODE profile is obtained by integrating
    \f[
        \frac{du^+}{dy^+} =
        \frac{1}{1 + \kappa y^+ \left(1 - \exp(-y^+/A^+)\right)^2}.
    \f]

    Evaluating the corresponding wall stress from a sampled velocity requires
    an iterative solve. This class replaces that solve by the explicit
    approximation of Nuca, Mukha, and Parsani. The approximation is constructed
    by adding a Gaussian perturbation to the explicit law of Cai and Sagaut.

    The approximation coefficients depend on \f$\kappa\f$ and \f$A^+\f$.
    The \c approximant entry controls which coefficient set is used:
    \verbatim
        auto       Select highRe for kappa=0.387, Aplus=15.2516;
                   select classical for kappa=0.41, Aplus=17;
                   otherwise select global.
        highRe     Fixed three-Gaussian fit for high-Re constants
                   kappa=0.387, Aplus=15.2516. The reported maximum relative
                   error with respect to the original equilibrium ODE model is
                   below 0.04%.
        classical  Fixed three-Gaussian fit for classical constants
                   kappa=0.41, Aplus=17. The reported maximum relative
                   error with respect to the original equilibrium ODE model is
                   below 0.04%.
        global     One-Gaussian regression in kappa and Aplus for other
                   constants. This is the fallback for user-selected constants
                   and is constructed to remain within the 1% relative-error
                   target.
    \endverbatim

    Usage:
    \verbatim
    Law
    {
        type        EquilibriumODE;
        kappa      value; (default 0.41)
        Aplus      value; (default 17)
        approximant auto; (default auto; highRe, classical or global)
    }
    \endverbatim

    Reference for the explicit approximation:
    \verbatim
        Nuca, R., Mukha, T., and Parsani, M. (2025).
        Explicit formulations of widely used wall models for large-eddy
        simulation.
        Physics of Fluids, 37, 035215.
        https://doi.org/10.1063/5.0253882
    \endverbatim

Contributors/Copyright:
    2024-2026 Timofey Mukha


SourceFiles
    EquilibriumODEExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef EquilibriumODEExplicitLawOfTheWall_H
#define EquilibriumODEExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"

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

        //- Equivalent log-law B constant for the ODE profile
        scalar B_;

        //- Approximation family: auto, highRe, classical or global
        word approximant_;

        //- Cai-Sagaut blending parameters
        scalar p_;
        scalar s_;

        //- Gaussian perturbation coefficients
        label nGaussians_;
        scalar mu_[3];
        scalar sigma_[3];
        scalar xi_[3];

        //- Select and store the approximation coefficients
        void setApproximantCoeffs(const word& approximant);

        //- Return whether two constants should be considered matching
        static bool approxEqual(const scalar a, const scalar b);

        //- Compute the Cai-Sagaut component of u+
        scalar CaiSagautUPlus(const scalar Re) const;

        //- Compute the perturbation component of u+
        scalar deltaUPlus(const scalar log10Re) const;

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

        //- Return the equivalent log-law B constant
        scalar B() const
        {
            return B_;
        }

        //- Return the selected approximant
        word approximant() const
        {
            return approximant_;
        }

        //- Return the Cai-Sagaut p coefficient
        scalar p() const
        {
            return p_;
        }

        //- Return the Cai-Sagaut s coefficient
        scalar s() const
        {
            return s_;
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
