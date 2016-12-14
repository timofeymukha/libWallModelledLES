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
    Newton root finder for better stability.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Aleksandar Jemcov.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "NewtonRoot.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::NewtonRoot::maxIter = 60;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRoot::NewtonRoot
(
    scalar (*f) (scalar),
    scalar (*d) (scalar),
    const scalar eps
)
:
    f_(f),
    d_(d),
    eps_(eps)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRoot::root
(
    scalar xOld
) const
{
    if (0 == d_(xOld))
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::NewtonRoot<Func, Deriv>::root\n"
            "(\n"
            "    const scalar xOld,\n"
            ") const"
        )   << "Jacobian equal to zero.  f'(xN) = " << d_(xOld)
            << abort(FatalError); }


    for (label nIter = 0; nIter < maxIter; ++nIter)
    {
        scalar f = this->f_(xOld);
        scalar d = this->d_(xOld);

        scalar xNew = xOld - f/d;

        if (mag(xNew - xOld) <= this->eps_)
        {
            return xNew;
        }

        xOld = xNew;
    }

    return xOld;
}


// ************************************************************************* //
