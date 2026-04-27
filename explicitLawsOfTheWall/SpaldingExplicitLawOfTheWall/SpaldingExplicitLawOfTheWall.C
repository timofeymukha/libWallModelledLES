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
    approximant_("auto"),
    p_(0),
    s_(0),
    nGaussians_(0),
    mu_{0, 0, 0},
    sigma_{0, 0, 0},
    xi_{0, 0, 0}
{
    constDict_.add("kappa", kappa);
    constDict_.add("B", B);
    setApproximantCoeffs(approximant_);
    constDict_.add("approximant", approximant_);

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

bool Foam::SpaldingExplicitLawOfTheWall::approxEqual
(
    const scalar a,
    const scalar b
)
{
    return mag(a - b) < 1e-4;
}


void Foam::SpaldingExplicitLawOfTheWall::setApproximantCoeffs
(
    const word& approximant
)
{
    word selected(approximant);

    if (selected == "auto")
    {
        if (approxEqual(kappa_, 0.387) && approxEqual(B_, 4.21))
        {
            selected = "highRe";
        }
        else if (approxEqual(kappa_, 0.4) && approxEqual(B_, 5.5))
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
        mu_[0] = 3.083144724251491; sigma_[0] = 2.3; xi_[0] = -0.09201260181506102;
        mu_[1] = 2.362020355550757; sigma_[1] = 2.298019985567246; xi_[1] = -0.32592545933705347;
        mu_[2] = 2.984943354098199; sigma_[2] = 0.5678593976861175; xi_[2] = 0.024978818026146963;
        p_ = 1.1768259148434161;
        s_ = 200.00359862759768;
    }
    else if (selected == "classical")
    {
        nGaussians_ = 3;
        mu_[0] = 3.072995395381741; sigma_[0] = 2.886186935612637; xi_[0] = 0.2888114606063557;
        mu_[1] = 2.5767807695580207; sigma_[1] = 0.9654690496618444; xi_[1] = -0.14752950593696862;
        mu_[2] = 2.623363478729046; sigma_[2] = 1.8670446391015763; xi_[2] = -1.2902885484930124;
        p_ = 1.1419065322343513;
        s_ = 400.0;
    }
    else if (selected == "global")
    {
        nGaussians_ = 1;
        mu_[0] = -1.8082*kappa_ + 0.0618*B_ + 2.4191;
        sigma_[0] = 0.3866*kappa_ + 0.0666*B_ + 0.7957;
        xi_[0] = -0.3558*kappa_ + 0.0386*B_ + 0.1890;
        mu_[1] = 0; sigma_[1] = 0; xi_[1] = 0;
        mu_[2] = 0; sigma_[2] = 0; xi_[2] = 0;
        p_ = 0.1161*kappa_ + 0.0203*B_ + 1.1276;
        s_ = -994.2937*kappa_ + 80.3579*B_ + 301.0408;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown Spalding explicit approximant " << selected << nl
            << "Valid options are auto, highRe, classical and global."
            << abort(FatalError);
    }

    approximant_ = selected;
}


Foam::scalar Foam::SpaldingExplicitLawOfTheWall::CaiSagautUPlus
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


Foam::scalar Foam::SpaldingExplicitLawOfTheWall::deltaUPlus
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


void Foam::SpaldingExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Explicit Spalding law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< indent << "approximant" << indent <<  approximant_ << nl;
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

    const scalar uPlus = CaiSagautUPlus(re) + deltaUPlus(Foam::log10(re));

    return u / uPlus;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
