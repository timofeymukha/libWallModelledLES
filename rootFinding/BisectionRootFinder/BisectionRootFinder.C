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
    scalar guess
) const
{
    label i = 1;
   
    scalar a = 1/bracket_*guess;
    scalar b  = bracket_*guess;
    scalar c = 0;
    scalar fC = 0;
    
    if (f_(a)*f_(b) >= 0)
    {
        // Increase interval range towards the wall
        a = SMALL;
    }

    if (f_(a)*f_(b) >= 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::BisectionRootFinder::root\n"
            "(\n"
            "  scalar guess\n"
            ") const"
        )   << "Root is not bracketed.  f(a) = " << f_(a) << " f(b) = " << f_(b)
            << abort(FatalError);
    }
    
    while (i <= maxIter_)
    {
       c = 0.5*(a + b);
       fC = f_(c);
       
       if ( (fC < SMALL) && (0.5*(b - a) < eps_) )
       {
           return c;
       }
       i++;
       
       if (sign(fC) == sign(f_(a)) )
       {
           a = c;          
       }
       else
       {
           b = c;     
       }
    }
   
    if (debug)
    {
        WarningIn
        (
            "Foam::scalar Foam::BisectionRootFinder::root\n"
            "(\n"
            "    scalar guess,\n"
            ") const"
        )   << "Maximum number of iterations exceeded";
    }
   
   return c;
}


// ************************************************************************* //
