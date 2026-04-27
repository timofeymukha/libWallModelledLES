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

    Reichardt's law is explicit in \f$u^+(y^+)\f$ but implicit in the wall
    stress when the sampled velocity and wall-normal distance are known. This
    class replaces the resulting nonlinear wall-stress solve by the explicit
    approximation of Nuca, Mukha, and Parsani. The approximation is constructed
    by adding a Gaussian perturbation to the explicit law of Cai and Sagaut.

    The approximation coefficients depend on the model constants. The
    \c approximant entry controls which coefficient set is used:
    \verbatim
        auto       Select highRe for kappa=0.387, B1=11, B2=3 and
                   C=6.663;
                   select classical for kappa=0.41, B1=11, B2=3 and C=7.8;
                   otherwise select global.
        highRe     Fixed three-Gaussian fit for the high-Re constants above.
                   The reported maximum relative error with respect to the
                   original Reichardt law is below 0.04%.
        classical  Fixed three-Gaussian fit for the classical constants above.
                   The reported maximum relative error with respect to the
                   original Reichardt law is below 0.04%.
        global     One-Gaussian regression in kappa and C for B1=11, B2=3.
                   This is the fallback for user-selected constants and is
                   constructed to remain within the 1% relative-error target.
    \endverbatim

    Usage:
    \verbatim
    Law
    {
        type      Reichardt;
        kappa     value; (default 0.41)
        B1        value; (default 11)
        B2        value; (default 3)
        C         value; (default 7.8)
        approximant auto; (default auto; highRe, classical or global)
    }
    \endverbatim

    Reference for the original law:
    \verbatim
        Reichardt, H. (1951).
        Vollstandige Darstellung der turbulenten Geschwindigkeitsverteilung in
        glatten Leitungen.
        Zeitschrift fur Angewandte Mathematik und Mechanik 31(7) (pp. 208-219).
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
    ReichardtExplicitLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef ReichardtExplicitLawOfTheWall_H
#define ReichardtExplicitLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "ExplicitLawOfTheWall.H"

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

        //- Equivalent logarithmic intercept, C + log(kappa)/kappa
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
