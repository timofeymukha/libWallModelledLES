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
    B_(C_ + Foam::log(kappa_)/kappa_),
    approximant_("auto"),
    p_(0),
    s_(0),
    nGaussians_(0),
    mu_{0, 0, 0},
    sigma_{0, 0, 0},
    xi_{0, 0, 0}
{
    constDict_.add("kappa", kappa);
    constDict_.add("B1", B1);
    constDict_.add("B2", B2);
    constDict_.add("C", C);
    setApproximantCoeffs(approximant_);
    constDict_.add("approximant", approximant_);

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
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.41)),
    B1_(constDict_.lookupOrAddDefault<scalar>("B1", 11)),
    B2_(constDict_.lookupOrAddDefault<scalar>("B2", 3)),
    C_(constDict_.lookupOrAddDefault<scalar>("C", 7.8)),
    B_(C_ + Foam::log(kappa_)/kappa_),
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

bool Foam::ReichardtExplicitLawOfTheWall::approxEqual
(
    const scalar a,
    const scalar b
)
{
    return mag(a - b) < 1e-4;
}


void Foam::ReichardtExplicitLawOfTheWall::setApproximantCoeffs
(
    const word& approximant
)
{
    word selected(approximant);

    if (selected == "auto")
    {
        if
        (
            approxEqual(kappa_, 0.387)
         && approxEqual(B1_, 11)
         && approxEqual(B2_, 3)
         && approxEqual(C_, 6.663)
        )
        {
            selected = "highRe";
        }
        else if
        (
            approxEqual(kappa_, 0.41)
         && approxEqual(B1_, 11)
         && approxEqual(B2_, 3)
         && approxEqual(C_, 7.8)
        )
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
        mu_[0] = 1.0; sigma_[0] = 0.6216482021221031; xi_[0] = -0.0024249886056448367;
        mu_[1] = 2.367258578988584; sigma_[1] = 2.4404184364171475; xi_[1] = 0.1213581213427751;
        mu_[2] = 2.397399420039065; sigma_[2] = 0.8459981526692487; xi_[2] = 0.08281299107063853;
        p_ = 1.2490442812025897;
        s_ = 131.58255333367103;
    }
    else if (selected == "classical")
    {
        nGaussians_ = 3;
        mu_[0] = 2.2998912899998016; sigma_[0] = 0.7189577648274934; xi_[0] = 0.09376484880628698;
        mu_[1] = 2.050645207323347; sigma_[1] = 1.7512288748770573; xi_[1] = -0.1;
        mu_[2] = 3.690238876998926; sigma_[2] = 1.5791921478874458; xi_[2] = -0.01401387317553293;
        p_ = 1.2331111626598046;
        s_ = 130.88152902954403;
    }
    else if (selected == "global")
    {
        nGaussians_ = 1;
        mu_[0] = 4.904772719171782*kappa_ - 0.3617181238520316*C_ + 2.907689416301654;
        sigma_[0] = 1.3041022247754417*kappa_ - 0.05465572961345442*C_ + 0.6762546999159017;
        xi_[0] = -1.0228920326709718*kappa_ + 0.04518022065539677*C_ + 0.16875743912554192;
        mu_[1] = 0; sigma_[1] = 0; xi_[1] = 0;
        mu_[2] = 0; sigma_[2] = 0; xi_[2] = 0;
        p_ = 0.3472185702403766*kappa_ - 0.004391948001728183*C_ + 1.1436513389571419;
        s_ = -99.55156119159305*kappa_ + 19.647991075591285*C_ + 18.173005459510563;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown Reichardt explicit approximant " << selected << nl
            << "Valid options are auto, highRe, classical and global."
            << abort(FatalError);
    }

    approximant_ = selected;
}


Foam::scalar Foam::ReichardtExplicitLawOfTheWall::CaiSagautUPlus
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


Foam::scalar Foam::ReichardtExplicitLawOfTheWall::deltaUPlus
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


void Foam::ReichardtExplicitLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Reichardt law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B1" << indent <<  B1_ << nl;
    Info<< indent << "B2" << indent <<  B2_ << nl;
    Info<< indent << "C" << indent <<  C_ << nl;
    Info<< indent << "approximant" << indent <<  approximant_ << nl;
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

    const scalar uPlus = CaiSagautUPlus(re) + deltaUPlus(Foam::log10(re));

    return u / uPlus;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
