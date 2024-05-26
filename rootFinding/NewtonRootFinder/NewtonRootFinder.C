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

#include "NewtonRootFinder.H"
#include "addToRunTimeSelectionTable.H"

using namespace boost::math::policies;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(NewtonRootFinder, 0);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, Word);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, Dictionary);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, DictionaryOnly);
}
#endif

typedef policy<
      domain_error<ignore_error>,
      overflow_error<ignore_error>
      > myPolicy;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRootFinder::root(
    scalar guess,
    scalar lowerBound,
    scalar upperBound
) const
{

    auto funcTuple = [this](scalar uTau)
    {
        return std::make_tuple(this->f_(uTau), this->d_(uTau));
    };

    auto maxIter = static_cast<boost::uintmax_t>(maxIter_);

    scalar result = boost::math::tools::newton_raphson_iterate(
        funcTuple,
        guess,
        lowerBound,
        upperBound,
        getDigits_,
        maxIter
    );

    if (debug)
    {
        WarningIn("Foam::NewtonRootFinder::root()")
            << "The method did not converge to desired tolerance." << nl;
    }

    return result;
}

// ************************************************************************* //
