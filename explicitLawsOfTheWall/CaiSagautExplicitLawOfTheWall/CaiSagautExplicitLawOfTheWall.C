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

\*---------------------------------------------------------------------------*/

#include "CaiSagautExplicitLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"
#include <boost/math/special_functions/lambert_w.hpp>

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CaiSagautExplicitLawOfTheWall, 0);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, CaiSagautExplicitLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, CaiSagautExplicitLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CaiSagautExplicitLawOfTheWall::CaiSagautExplicitLawOfTheWall
(
    const scalar kappa,
    const scalar B,
    const scalar p,
    const scalar s
)
:
    ExplicitLawOfTheWall(),
    kappa_(kappa),
    B_(B),
    p_(p),
    s_(s)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B", B);
    constDict_.add("p", p);
    constDict_.add("s", s);

    if (debug)
    {
        printCoeffs();
    }

}

Foam::CaiSagautExplicitLawOfTheWall::CaiSagautExplicitLawOfTheWall
(
    const dictionary & dict
)
:
    ExplicitLawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    B_(constDict_.lookupOrAddDefault<scalar>("B", 5.5)),
    p_(constDict_.lookupOrAddDefault<scalar>("p", 1.138)),
    s_(constDict_.lookupOrAddDefault<scalar>("s", 217.8))
{
    if (debug)
    {
        printCoeffs();
    }

}

Foam::CaiSagautExplicitLawOfTheWall::CaiSagautExplicitLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    CaiSagautExplicitLawOfTheWall(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CaiSagautExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "CaiSagaut law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< indent << "p" << indent <<  p_ << nl;
    Info<< indent << "s" << indent <<  s_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::CaiSagautExplicitLawOfTheWall::uTau
(
    const SingleCellSampler & sampler,
    label index,
    scalar nu
) const
{
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");
    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    scalar y = sampler.h()[index];

    const scalar re = u * y / nu;
    const scalar f = exp(-re / s_);
    const scalar E = exp(kappa_ * B_);

    scalar uPlus = Foam::pow(f, p_) * Foam::sqrt(re) ;
    uPlus += Foam::pow(1 - f, p_) / kappa_ * boost::math::lambert_w0(kappa_*E*re);

    return u / uPlus;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
