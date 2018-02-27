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

#include "ReichardtLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(ReichardtLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, ReichardtLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(LawOfTheWall, ReichardtLawOfTheWall, TypeAndDictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ReichardtLawOfTheWall::ReichardtLawOfTheWall
(
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4)),
    B1_(dict.lookupOrDefault<scalar>("B1", 11)),
    B2_(dict.lookupOrDefault<scalar>("B2", 3)),
    C_(dict.lookupOrDefault<scalar>("C", 7.8))
{
    if (debug)
    {        
        printCoeffs();
    }
    
}

Foam::ReichardtLawOfTheWall::ReichardtLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4)),
    B1_(dict.lookupOrDefault<scalar>("B1", 11)),
    B2_(dict.lookupOrDefault<scalar>("B2", 3)),
    C_(dict.lookupOrDefault<scalar>("C", 7.8))
{
    if (debug)
    {        
        printCoeffs();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ReichardtLawOfTheWall::printCoeffs() const
{
    Info<< nl << "Reichardt law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B1" << indent <<  B1_ << nl;
    Info<< indent << "B2" << indent <<  B2_ << nl;
    Info<< indent << "C" << indent <<  C_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::ReichardtLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{  
    scalar uPlus = u/uTau;
    scalar yPlus = y*uTau/nu;

    return uPlus - 1/kappa_*log(1 + kappa_*yPlus) -
           C_*(1 - exp(-yPlus/B1_) -
           yPlus/B1_*exp(-yPlus/B2_));
}


Foam::scalar Foam::ReichardtLawOfTheWall::derivative
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu        
) const
{
    scalar uPlus = u/uTau;
    scalar yPlus = y*uTau/nu;
    
    return -uPlus/uTau - y/nu/(1 + kappa_*yPlus) - 
           C_*y/(nu*B1_)*(exp(-yPlus/B1_) -
           exp(-yPlus/B2_) + yPlus/B2_*exp(-yPlus/B2_));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //