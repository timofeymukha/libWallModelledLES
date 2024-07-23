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

#include "ReichardtExplicitLawOfTheWall.H"
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
    defineTypeNameAndDebug(ReichardtExplicitLawOfTheWall, 0);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, ReichardtExplicitLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, ReichardtExplicitLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ReichardtExplicitLawOfTheWall::ReichardtExplicitLawOfTheWall
(
    const scalar kappa,
    const scalar B1,
    const scalar B2,
    const scalar C
)
:
    ExplicitLawOfTheWall(),
    kappa_(kappa),
    B1_(B1),
    B2_(B2),
    C_(C),
    CaiSagaut_(kappa, C + Foam::log(kappa)/kappa, 1.24115752, 117.36295084)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B1", B1);
    constDict_.add("B2", B2);
    constDict_.add("C", C);

    if (debug)
    {
        printCoeffs();
    }

}

Foam::ReichardtExplicitLawOfTheWall::ReichardtExplicitLawOfTheWall
(
    const dictionary & dict
)
:
    ExplicitLawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    B1_(constDict_.lookupOrAddDefault<scalar>("B1", 11)),
    B2_(constDict_.lookupOrAddDefault<scalar>("B2", 3)),
    C_(constDict_.lookupOrAddDefault<scalar>("C", 7.8)),
    CaiSagaut_(kappa_, C_ + Foam::log(kappa_)/kappa_, 1.24115752, 117.36295084)

{
    if (debug)
    {
        printCoeffs();
    }
}

Foam::ReichardtExplicitLawOfTheWall::ReichardtExplicitLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    ReichardtExplicitLawOfTheWall(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ReichardtExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Reichardt law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B1" << indent <<  B1_ << nl;
    Info<< indent << "B2" << indent <<  B2_ << nl;
    Info<< indent << "C" << indent <<  C_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::ReichardtExplicitLawOfTheWall::uTau
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

    scalar uPlus = CaiSagaut_.uTau(sampler, index, nu);
    uPlus += Helpers::gaussian(2.21217354, 0.51791062, -0.17108678, Foam::log10(re));
    uPlus += Helpers::gaussian(2.16739418, 0.88976565, 0.22467868, Foam::log10(re));
    uPlus += Helpers::gaussian(2.99178237, 0.20894151, 0.02020086, Foam::log10(re));

    return u / uPlus;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
