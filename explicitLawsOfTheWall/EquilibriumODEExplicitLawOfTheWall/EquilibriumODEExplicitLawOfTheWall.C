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

#include "EquilibriumODEExplicitLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"
#include <boost/math/special_functions/lambert_w.hpp>
#include "helpers.H"
#include "AdaptiveIntegrator.hpp"
#include <functional>

typedef std::function<Foam::scalar(const Foam::scalar)> IntegrandFunc;

Foam::scalar get_B(const Foam::scalar kappa, const Foam::scalar Aplus);

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(EquilibriumODEExplicitLawOfTheWall, 0);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, EquilibriumODEExplicitLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(ExplicitLawOfTheWall, EquilibriumODEExplicitLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EquilibriumODEExplicitLawOfTheWall::EquilibriumODEExplicitLawOfTheWall
(
    const scalar kappa,
    const scalar Aplus
)
:
    ExplicitLawOfTheWall(),
    kappa_(kappa),
    Aplus_(Aplus),
    CaiSagaut_(kappa_, 0, 1.19072560e+00, 2.18729438e+02)
{
    constDict_.add("kappa", kappa);
    constDict_.add("Aplus", Aplus);

//            nutp = lambda yp, kappa, Aplus: kappa*yp*(1 - np.exp(-yp/Aplus))**2
//        integrand = lambda yp, kappa, Aplus: 1/(1 + nutp(yp, kappa, Aplus))
//        return np.array([quad(integrand, 0, val, args=(kappa, Aplus))[0] for val in yp])

    const scalar B =  get_B(kappa, Aplus);
    this->CaiSagaut_.set_B(4.21);

    Info << "B" << B << nl;

    if (debug)
    {
        printCoeffs();
    }

}

Foam::EquilibriumODEExplicitLawOfTheWall::EquilibriumODEExplicitLawOfTheWall
(
    const dictionary & dict
)
:
    ExplicitLawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.41)),
    Aplus_(constDict_.lookupOrAddDefault<scalar>("Aplus", 17)),
    CaiSagaut_(kappa_, 0, 1.19072560e+00, 2.18729438e+02)

{
    const scalar B =  get_B(this->kappa(), this->Aplus());
    this->CaiSagaut_.set_B(4.21);

    Info << "!!!!!!!!!!!!!!!!!!!!!!! B" << B << nl;

    if (debug)
    {
        printCoeffs();
    }
}

Foam::EquilibriumODEExplicitLawOfTheWall::EquilibriumODEExplicitLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    EquilibriumODEExplicitLawOfTheWall(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EquilibriumODEExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "EquilibriumODEExplicit law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "Aplus" << indent <<  Aplus_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::EquilibriumODEExplicitLawOfTheWall::uTau
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
    scalar uPlusDelta = 0;

    uPlusDelta += Helpers::gaussian(2.74758381e+00, 2.98641869e+00,  1.41579957e-01, Foam::log10(re));
    uPlusDelta += Helpers::gaussian(4.15537857e+00, 9.45524229e-01, -3.12760808e-02, Foam::log10(re));
    uPlusDelta += Helpers::gaussian(3.21739839e+00, 5.81214845e-01,  5.31676710e-02, Foam::log10(re));

    return u / (uPlus + uPlusDelta);

}

Foam::scalar get_B(const Foam::scalar kappa, const Foam::scalar Aplus)
{
    IntegrandFunc nut_func = [kappa, Aplus](const Foam::scalar yPlus)
    {
        return kappa * yPlus * Foam::sqr(1 - Foam::exp(-yPlus / Aplus));
    };

    IntegrandFunc integrand = [nut_func](const Foam::scalar yPlus)
    {
        return 1.0 / (1.0 + nut_func(yPlus));
    };

    Foam::scalar yPlus = 100000;

    AdaptiveIntegrator<Foam::scalar (Foam::scalar)> quad;
    Foam::scalar up = quad.integrate(integrand, 0.0, yPlus, 1e-3);

    return up - 1/kappa * Foam::log(yPlus);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
