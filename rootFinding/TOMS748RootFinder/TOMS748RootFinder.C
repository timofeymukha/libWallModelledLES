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

#include "TOMS748RootFinder.H"
#include "addToRunTimeSelectionTable.H"

using namespace boost::math::policies;
using namespace boost::math::tools;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(TOMS748RootFinder, 0);
    addToRunTimeSelectionTable(RootFinder, TOMS748RootFinder, Word);
    addToRunTimeSelectionTable(RootFinder, TOMS748RootFinder, Dictionary);
    addToRunTimeSelectionTable(RootFinder, TOMS748RootFinder, DictionaryOnly);
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::pair<Foam::scalar, Foam::label> Foam::TOMS748RootFinder::root(
    scalar guess,
    scalar lowerBound,
    scalar upperBound
) const
{

    auto maxIter = static_cast<boost::uintmax_t>(maxIter_);
    eps_tolerance<scalar> tol(getDigits_);

    label iterations = 0;
    auto wrapper = [this, &iterations](scalar uTau)
    {
        iterations++;
        return f_(uTau);
    };

    std::pair<scalar, scalar> result =
        toms748_solve(wrapper, lowerBound, upperBound, tol, maxIter);

    return std::make_pair(0.5*(result.first + result.second), iterations);
}

// ************************************************************************* //
