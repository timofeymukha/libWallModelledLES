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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRootFinder::root(scalar guess) const
{
    scalar error = 0;

    // Crash if derivative is zero
    if (0 == d_(guess))
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::NewtonRoot::root\n"
            "(\n"
            "    scalar guess,\n"
            ") const"
        )   << "Derivative equal to zero.  f' = " << d_(guess)
            << abort(Foam::FatalError); }

    for (label nIter = 0; nIter < maxIter_; ++nIter)
    {
        scalar f = f_(guess);
        scalar d = d_(guess);

        scalar newGuess = guess - f/d;
        
        error = mag(newGuess - guess)/mag(guess);
    
        if (error <= eps_)
        {
            return newGuess;
        }
        
        guess = newGuess;
    }
    
    if (debug)
    {
        WarningIn("Foam::NewtonRootFinder::root()")
            << "The method did not converge to desired tolerance." << nl;
    }
    return guess;
}

// ************************************************************************* //
