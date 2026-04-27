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

#include "BisectionRootFinder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(BisectionRootFinder, 0);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, Word);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, Dictionary);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, DictionaryOnly);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BisectionRootFinder::BisectionRootFinder
(
    const word & rootFinderName,
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const label maxIter
)
:
    RootFinder(rootFinderName, f, d, maxIter)
{}

Foam::BisectionRootFinder::BisectionRootFinder
(
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary & dict
)
:
    RootFinder(f, d, dict)
{}

Foam::BisectionRootFinder::BisectionRootFinder
(
            const dictionary & dict
)
:
    RootFinder(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::pair<Foam::scalar, Foam::label> Foam::BisectionRootFinder::root
(
    scalar guess,
    scalar lowerBound,
    scalar upperBound
) const
{
    scalar a = lowerBound;
    scalar b = upperBound;
    scalar fA = f_(a);
    scalar fB = f_(b);

    if (a > b)
    {
        Foam::Swap(a, b);
        Foam::Swap(fA, fB);
    }

    if (mag(fA) <= ROOTVSMALL)
    {
        return std::make_pair(a, 0);
    }

    if (mag(fB) <= ROOTVSMALL)
    {
        return std::make_pair(b, 0);
    }

    if (sign(fA) == sign(fB))
    {
        FatalErrorIn
        (
            "Foam::BisectionRootFinder::root\n"
            "(\n"
            "    scalar guess,\n"
            "    scalar lowerBound,\n"
            "    scalar upperBound\n"
            ") const"
        )   << "Root is not bracketed. f(" << a << ") = " << fA
            << " and f(" << b << ") = " << fB
            << abort(FatalError);
    }

    boost::uintmax_t maxIter = static_cast<boost::uintmax_t>(maxIter_);
    boost::math::tools::eps_tolerance<scalar> tolerance(getDigits_);

    label evaluations = 0;
    auto wrapper = [this, &evaluations](const scalar x)
    {
        evaluations++;
        return f_(x);
    };

    std::pair<scalar, scalar> result =
        boost::math::tools::bisect(wrapper, a, b, tolerance, maxIter);

    if (debug)
    {
        if (maxIter >= static_cast<boost::uintmax_t>(maxIter_))
        {
            WarningIn
            (
                "Foam::BisectionRootFinder::root\n"
                "(\n"
                "    scalar guess,\n"
                "    scalar lowerBound,\n"
                "    scalar upperBound\n"
                ") const"
            )   << "Maximum number of iterations exceeded";
        }
    }

    return std::make_pair(0.5*(result.first + result.second), evaluations);
}


// ************************************************************************* //
