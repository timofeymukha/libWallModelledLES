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

#include "NewtonRoot.H"
#include "error.H"
#include <functional>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRoot::NewtonRoot
(
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const scalar eps,
    const label maxIter
)
:
    f_(f),
    d_(d),
    eps_(eps),
    maxIter_(maxIter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRoot::root
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


    //Info<< "Newton initial" << guess << endl;
    for (label nIter = 0; nIter < maxIter_; ++nIter)
    {
        scalar f = this->f_(guess);
        scalar d = this->d_(guess);

        scalar xNew = guess - f/d;

     /*   if (mag(xNew - guess) <= this->eps_)
        {
      //      Info<< "Newton guess" << guess << " " << nIter << endl;
      //      Info<< "Newton new" << xNew << " " << nIter << endl;
            return xNew;
        } */
        
        guess = xNew;
       // Info<< guess << endl;
    }

    return guess;
}


// ************************************************************************* //
