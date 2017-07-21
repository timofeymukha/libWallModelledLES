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

namespace Foam
{
    defineTypeNameAndDebug(BisectionRootFinder, 0);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, Word);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, Dictionary);
    addToRunTimeSelectionTable(RootFinder, BisectionRootFinder, DictionaryOnly);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BisectionRootFinder::BisectionRootFinder
(
    const word & rootFinderName, 
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const scalar eps,
    const label maxIter
)
:
    RootFinder(rootFinderName, f, d, eps, maxIter),
    bracket_(3.0)
{}

Foam::BisectionRootFinder::BisectionRootFinder
(
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary & dict
)
:
    RootFinder(f, d, dict),
    bracket_(dict.lookupOrDefault<scalar>("bracket", 3.0))
{}

Foam::BisectionRootFinder::BisectionRootFinder
(
            const dictionary & dict
)
:
    RootFinder(dict),
    bracket_(dict.lookupOrDefault<scalar>("bracket", 3.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// This one probably needs some work
Foam::scalar Foam::BisectionRootFinder::root
(
    scalar x1
) const
{
    scalar f, fMid, dx, rtb, xMid;

    // Create interval based on bracket and initial guess
    scalar x0 = 1/bracket_*x1;
    x1 = bracket_*x1;
    
    fMid = f_(x1);  
    f = f_(x0);
    
  
    // Crash if root is not bracketed. We might want to fall back to
    // [0, bracket*x1] first, since uTau > 0. This may be needed in the first
    // iterations of the solver
    if (f*fMid >= 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::BisectionRootFinder::root\n"
            "(\n"
            "  scalar x1\n"
            ") const"
        )   << "Root is not bracketed.  f(x0) = " << f << " f(x1) = " << fMid
            << abort(FatalError);
    }

    // Orient the search such that f > 0 lies at x + dx
    if (f < 0)
    {
        dx = x1 - x0;
        rtb = x0;
    }
    else
    {
        dx = x0 - x1;
        rtb = x1;
    }

    for (label nIter = 0; nIter < maxIter_; nIter++)
    {
        dx *= 0.5;
        xMid = rtb + dx;

        fMid = f_(xMid);

        if (fMid <= 0)
        {
            rtb = xMid;
        }

        // Note that this is the absolute error, and we use relative in Newton..
        if (mag(dx)/rtb < eps_ || mag(fMid) < SMALL)
        {
            return rtb;
        }
    }

    if (debug)
    {
        WarningIn
        (
            "Foam::scalar Foam::BisectionRoot<Func>::root\n"
            "(\n"
            "    const scalar x0,\n"
            "    const scalar x1\n"
            ") const"
        )   << "Maximum number of iterations exceeded";
    }
    return rtb;
}


// ************************************************************************* //
