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

Class
    NewtonRoot

Description
    Newton root finder.

Author
    Timofey Mukha.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "NewtonRootFinder.H"
#include "RootFinder.H"
#include "word.H"
#include "error.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include <functional>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NewtonRootFinder, 0);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, Word);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, Dictionary);
    addToRunTimeSelectionTable(RootFinder, NewtonRootFinder, DictionaryOnly);
}

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
            "    scalar x0,\n"
            ") const"
        )   << "Jacobian equal to zero.  f'(xN) = " << d_(guess)
            << abort(FatalError); }

    for (label nIter = 0; nIter < maxIter_; ++nIter)
    {
        scalar f = f_(guess);
        scalar d = d_(guess);

        scalar xNew = guess - f/d;
        
        error = mag(xNew - guess)/mag(guess);
    
        if (error <= eps_)
        {
            return xNew;
        }
        
        guess = xNew;
    }
    
    if (debug)
    {
        WarningIn("Foam::NewtonRootFinder::root()")
            << "The method did not converge to desired tolerance." << nl;
    }
    return guess;
}

// ************************************************************************* //
