/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    NewtonRoot

Description
    Newton root finder.

Author
    Timofey Mukha.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "NewtonRootFinder.H"
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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRootFinder::NewtonRootFinder
(
    const word & rootFinderName, 
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const scalar eps,
    const label maxIter
)
:
    RootFinder(rootFinderName, f, d, eps, maxIter),
    d_(d)
{}

Foam::NewtonRootFinder::NewtonRootFinder
(
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const dictionary & dict
)
:
    RootFinder(f, d, dict),
    d_(d)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRootFinder::root
(
    scalar x0
) const
{
    scalar guess = x0;
    if (0 == d_(guess))
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::NewtonRoot::root\n"
            "(\n"
            "    scalar xOld,\n"
            ") const"
        )   << "Jacobian equal to zero.  f'(xN) = " << d_(guess)
            << abort(FatalError); }

    for (label nIter = 0; nIter < maxIter_; ++nIter)
    {
        scalar f = this->f_(guess);
        scalar d = this->d_(guess);

        scalar xNew = guess - f/d;

        if (mag(xNew - guess)/mag(guess) <= this->eps_)
        {
            return xNew;
        }
        
        guess = xNew;
    }
    
    WarningIn("Foam::NewtonRootFinder::root()")
        << "The method did not converge to desired tolerance." << nl;

    return guess;
}


// ************************************************************************* //
