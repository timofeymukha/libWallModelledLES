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

#include "IntegratedReichardtLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(IntegratedReichardtLawOfTheWall, 0);
    addToRunTimeSelectionTable
    (
        LawOfTheWall,
        IntegratedReichardtLawOfTheWall,
        Dictionary
    );
    addToRunTimeSelectionTable
    (
        LawOfTheWall,
        IntegratedReichardtLawOfTheWall,
        TypeAndDictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IntegratedReichardtLawOfTheWall::IntegratedReichardtLawOfTheWall
(
    const dictionary & dict,
    const CellIndexList & list
)
:
    LawOfTheWall(dict, list),
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

Foam::IntegratedReichardtLawOfTheWall::IntegratedReichardtLawOfTheWall
(
    const word & lawName,
    const dictionary & dict,
    const CellIndexList & list
)
:
    LawOfTheWall(lawName, dict, list),
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

void Foam::IntegratedReichardtLawOfTheWall::printCoeffs() const
{
    Info<< nl << "IntegratedReichardt law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B1" << indent <<  B1_ << nl;
    Info<< indent << "B2" << indent <<  B2_ << nl;
    Info<< indent << "C" << indent <<  C_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::value
(
    scalar u,
    scalar index,
    scalar uTau,
    scalar nu
) const
{  
    scalar h = cellIndexList_.h()[index];
    
    scalar h1 = mag(h - cellIndexList_.lengthList()[index]/2);
    scalar h2 = h + cellIndexList_.lengthList()[index]/2;
    
    return u*(h2 - h1) - (logTerm(h2, uTau, nu) - logTerm(h1, uTau, nu) +
                          expTerm(h2, uTau, nu) - expTerm(h1, uTau, nu));
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::derivative
(
    scalar u,
    scalar index,
    scalar uTau,
    scalar nu        
) const
{
    scalar h = cellIndexList_.h()[index];
    
    scalar h1 = mag(h - cellIndexList_.lengthList()[index]/2);
    scalar h2 = h + cellIndexList_.lengthList()[index]/2;
        
    return -(logTermDerivative(h2, uTau, nu) - logTermDerivative(h1, uTau, nu) +
             expTermDerivative(h2, uTau, nu) - logTermDerivative(h1, uTau, nu));
}

Foam::scalar Foam::IntegratedReichardtLawOfTheWall::logTerm
(
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar yPlus = y*uTau/nu;
    
    return nu/kappa_*(-yPlus + log(1 + kappa_*yPlus)*(yPlus + 1/kappa_));
    
}

Foam::scalar Foam::IntegratedReichardtLawOfTheWall::expTerm
(
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar yPlus = y*uTau/nu;
    
    scalar term1 = yPlus;
    scalar term2 = B1_*exp(-yPlus/B1_);
    scalar term3 = B2_*(B2_ + yPlus)/B1_*exp(-yPlus/B2_);
    return C_*nu*(term1 + term2 + term3);
    
}

Foam::scalar Foam::IntegratedReichardtLawOfTheWall::logTermDerivative
(
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar yPlus = y*uTau/nu;
            
    return y*(yPlus/(kappa_*yPlus + 1) -
           1/kappa_ +
           1/kappa_*log(kappa_*yPlus + 1) +
           1/(kappa_*(kappa_*yPlus + 1)));    
}

Foam::scalar Foam::IntegratedReichardtLawOfTheWall::expTermDerivative
(
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar yPlus = y*uTau/nu;
    
    return C_*(y - y*exp(-yPlus/B1_) - y*yPlus/B1_*exp(-yPlus/B2_));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //