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

#include "CaiSagautLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"
#include <boost/math/special_functions/lambert_w.hpp>

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CaiSagautLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, CaiSagautLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(LawOfTheWall, CaiSagautLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CaiSagautLawOfTheWall::CaiSagautLawOfTheWall
(
    const scalar kappa,
    const scalar B
)
:
    LawOfTheWall(),
    kappa_(kappa),
    B_(B)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B", B);

    if (debug)
    {
        printCoeffs();
    }

}

Foam::CaiSagautLawOfTheWall::CaiSagautLawOfTheWall
(
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    B_(constDict_.lookupOrAddDefault<scalar>("B", 5.5))
{
    if (debug)
    {
        printCoeffs();
    }

}

Foam::CaiSagautLawOfTheWall::CaiSagautLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    CaiSagautLawOfTheWall(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CaiSagautLawOfTheWall::printCoeffs() const
{
    Info<< nl << "CaiSagaut law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::CaiSagautLawOfTheWall::value
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");
    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    scalar y = sampler.h()[index];

    return value(u, y, uTau, nu);
}

Foam::scalar Foam::CaiSagautLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    const scalar p = 1.138;
    const scalar s = 217.8;
    const scalar re = u * y / nu;
    const scalar f = exp(-re / s);
    const scalar E = exp(kappa_ * B_);

    scalar uPlus = Foam::pow(f, p) * Foam::sqrt(re) ;
    uPlus += Foam::pow(1 - f, p) / kappa_ * boost::math::lambert_w0(kappa_*E*re);

    return uTau - u / (uPlus + VSMALL);
}


Foam::scalar Foam::CaiSagautLawOfTheWall::derivative
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{
    return 1.0;
}


Foam::scalar Foam::CaiSagautLawOfTheWall::derivative
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{

    return 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
