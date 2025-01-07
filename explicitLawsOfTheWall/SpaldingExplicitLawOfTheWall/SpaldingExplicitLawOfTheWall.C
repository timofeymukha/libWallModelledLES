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

#include "SpaldingExplicitLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"
#include <boost/math/special_functions/lambert_w.hpp>
#include "helpers.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SpaldingExplicitLawOfTheWall, 0);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, SpaldingExplicitLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, SpaldingExplicitLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SpaldingExplicitLawOfTheWall::SpaldingExplicitLawOfTheWall
(
    const scalar kappa,
    const scalar B
)
:
    ExplicitLawOfTheWall(),
    kappa_(kappa),
    B_(B),
    CaiSagaut_(kappa_, B_, 1.1764, 206.0388)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B", B);

    if (debug)
    {
        printCoeffs();
    }

}
Foam::SpaldingExplicitLawOfTheWall::SpaldingExplicitLawOfTheWall
(
    const dictionary & dict
)
:
    ExplicitLawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    B_(constDict_.lookupOrAddDefault<scalar>("B", 5.5)),
    CaiSagaut_(kappa_, B_, 1.1764, 206.0388)

{
    if (debug)
    {
        printCoeffs();
    }
}

Foam::SpaldingExplicitLawOfTheWall::SpaldingExplicitLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    SpaldingExplicitLawOfTheWall(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SpaldingExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Explicit Spalding law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::SpaldingExplicitLawOfTheWall::uTau
(
    const SingleCellSampler & sampler,
    label index,
    scalar nu
) const
{
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");
    const scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    const scalar y = sampler.h()[index];
    const scalar re = u * y / nu;

    scalar uPlus = u / CaiSagaut_.uTau(sampler, index, nu);
    uPlus += Helpers::gaussian(3.1133, 2.3485, -0.0864, Foam::log10(re));
    uPlus += Helpers::gaussian(2.369, 2.274, -0.3484,  Foam::log10(re));
    uPlus += Helpers::gaussian(2.9937, 0.5698, 0.0254, Foam::log10(re));

    return u / uPlus;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
