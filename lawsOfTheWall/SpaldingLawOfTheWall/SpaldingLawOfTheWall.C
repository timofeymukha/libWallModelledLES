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
    SpaldingLawOfTheWall

Description
    The law of the wall proposed by Spalding.

Authors
    Timofey Mukha.  All rights reserved.

 \*---------------------------------------------------------------------------*/

#include "SpaldingLawOfTheWall.H"
#include "error.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SpaldingLawOfTheWall::SpaldingLawOfTheWall()
:
    kappa(0.4),
    B(5.5)
{}


Foam::SpaldingLawOfTheWall::SpaldingLawOfTheWall
(
    scalar kappa,
    scalar B
)
:
    kappa(kappa),
    B(B)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::SpaldingLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar uPlus = u/uTau;
    return uPlus + exp(-kappa*B)*(exp(kappa*uPlus) - 1 - kappa*uPlus
         - 0.5*sqr(kappa*uPlus) - 1./6*kappa*uPlus*sqr(kappa*uPlus))
         - y*uTau/nu;
}

Foam::scalar Foam::SpaldingLawOfTheWall::derivativeValue
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu        
) const
{
    scalar uPlus = u/uTau;
    return -y/nu - u/sqr(uTau) - kappa*uPlus/uTau*exp(-kappa*B)
           *(exp(kappa*uPlus) - 1 - kappa*uPlus - 0.5*sqr(kappa*uPlus));
}
// ************************************************************************* //
