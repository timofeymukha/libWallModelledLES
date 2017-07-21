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

#include "SpaldingLawOfTheWall.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(SpaldingLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, SpaldingLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(LawOfTheWall, SpaldingLawOfTheWall, TypeAndDictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SpaldingLawOfTheWall::SpaldingLawOfTheWall
(
    const Foam::dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4)),
    B_(dict.lookupOrDefault<scalar>("B", 5.5))
{
    if (debug)
    {        
        printCoeffs();
    }
}


Foam::SpaldingLawOfTheWall::SpaldingLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4)),
    B_(dict.lookupOrDefault<scalar>("B", 5.5))
{
    if (debug)
    {        
        printCoeffs();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SpaldingLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Spalding law of the wall" << nl;     
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}

Foam::scalar Foam::SpaldingLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar uPlus = u/uTau;
    return uPlus + exp(-kappa_*B_)*(exp(kappa_*uPlus) - 1 - kappa_*uPlus
         - 0.5*sqr(kappa_*uPlus) - 1./6*kappa_*uPlus*sqr(kappa_*uPlus))
         - y*uTau/nu;
}

Foam::scalar Foam::SpaldingLawOfTheWall::derivative
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu        
) const
{
    scalar uPlus = u/uTau;
    return -y/nu - u/sqr(uTau) - kappa_*uPlus/uTau*exp(-kappa_*B_)
           *(exp(kappa_*uPlus) - 1 - kappa_*uPlus - 0.5*sqr(kappa_*uPlus));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

