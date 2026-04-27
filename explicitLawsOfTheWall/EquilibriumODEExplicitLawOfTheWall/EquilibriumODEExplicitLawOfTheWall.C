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
    B_(get_B(kappa_, Aplus_)),
    approximant_("auto"),
    p_(0),
    s_(0),
    nGaussians_(0),
    mu_{0, 0, 0},
    sigma_{0, 0, 0},
    xi_{0, 0, 0}
{
    constDict_.add("kappa", kappa);
    constDict_.add("Aplus", Aplus);
    setApproximantCoeffs(approximant_);
    constDict_.add("approximant", approximant_);

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
    B_(get_B(kappa_, Aplus_)),
    approximant_(constDict_.lookupOrAddDefault<word>("approximant", "auto")),
    p_(0),
    s_(0),
    nGaussians_(0),
    mu_{0, 0, 0},
    sigma_{0, 0, 0},
    xi_{0, 0, 0}

{
    setApproximantCoeffs(approximant_);
    constDict_.set("approximant", approximant_);

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

bool Foam::EquilibriumODEExplicitLawOfTheWall::approxEqual
(
    const scalar a,
    const scalar b
)
{
    return mag(a - b) < 1e-6;
}


void Foam::EquilibriumODEExplicitLawOfTheWall::setApproximantCoeffs
(
    const word& approximant
)
{
    word selected(approximant);

    if (selected == "auto")
    {
        if (approxEqual(kappa_, 0.387) && approxEqual(Aplus_, 15.2516))
        {
            selected = "highRe";
        }
        else if (approxEqual(kappa_, 0.41) && approxEqual(Aplus_, 17))
        {
            selected = "classical";
        }
        else
        {
            selected = "global";
        }
    }

    if (selected == "highRe")
    {
        nGaussians_ = 3;
        mu_[0] = 2.747578886341323; sigma_[0] = 2.986331718406963; xi_[0] = 0.14158434859130883;
        mu_[1] = 3.217234942007338; sigma_[1] = 0.5812144329347989; xi_[1] = 0.053147124312226554;
        mu_[2] = 4.155551202580644; sigma_[2] = 0.945720164071122; xi_[2] = -0.03125681152478188;
        p_ = 1.1907235645597454;
        s_ = 218.72498935725267;
    }
    else if (selected == "classical")
    {
        nGaussians_ = 3;
        mu_[0] = 2.712195304468856; sigma_[0] = 1.153331759918931; xi_[0] = 0.04253708716912398;
        mu_[1] = 2.785538268044812; sigma_[1] = 3.1614868745529012; xi_[1] = 0.13731053348456193;
        mu_[2] = 2.547620704489551; sigma_[2] = 0.6124242434703746; xi_[2] = 0.012277142724028153;
        p_ = 1.2029211749614945;
        s_ = 247.0697345675632;
    }
    else if (selected == "global")
    {
        nGaussians_ = 1;
        mu_[0] = 3.498902867914008*kappa_ + 0.13030217331542102*Aplus_ - 0.9085784915858164;
        sigma_[0] = -1.2191000776979632*kappa_ - 0.02844761945101745*Aplus_ + 1.5494923375953826;
        xi_[0] = -0.39993239380818907*kappa_ - 0.0072855142653190505*Aplus_ + 0.3209799815846869;
        mu_[1] = 0; sigma_[1] = 0; xi_[1] = 0;
        mu_[2] = 0; sigma_[2] = 0; xi_[2] = 0;
        p_ = 0.22127889108172996*kappa_ + 0.0052185898765612845*Aplus_ + 1.0616974939891426;
        s_ = -130.48555166521916*kappa_ + 7.964453719974989*Aplus_ + 10.6631049793274;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown EquilibriumODE explicit approximant " << selected << nl
            << "Valid options are auto, highRe, classical and global."
            << abort(FatalError);
    }

    approximant_ = selected;
}


Foam::scalar Foam::EquilibriumODEExplicitLawOfTheWall::CaiSagautUPlus
(
    const scalar Re
) const
{
    const scalar f = exp(-Re / s_);
    const scalar E = exp(kappa_ * B_);

    scalar uPlus = Foam::pow(f, p_) * Foam::sqrt(Re);
    uPlus += Foam::pow(1 - f, p_) / kappa_
        * boost::math::lambert_w0(kappa_*E*Re);

    return uPlus;
}


Foam::scalar Foam::EquilibriumODEExplicitLawOfTheWall::deltaUPlus
(
    const scalar log10Re
) const
{
    scalar delta = 0;

    for (label i = 0; i < nGaussians_; i++)
    {
        delta += Helpers::gaussian(mu_[i], sigma_[i], xi_[i], log10Re);
    }

    return delta;
}


void Foam::EquilibriumODEExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "EquilibriumODEExplicit law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "Aplus" << indent <<  Aplus_ << nl;
    Info<< indent << "approximant" << indent <<  approximant_ << nl;
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

    const scalar uPlus = CaiSagautUPlus(re) + deltaUPlus(Foam::log10(re));

    return u / uPlus;

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
