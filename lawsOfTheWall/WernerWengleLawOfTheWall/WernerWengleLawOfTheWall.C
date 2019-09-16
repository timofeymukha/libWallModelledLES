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

#include "WernerWengleLawOfTheWall.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(WernerWengleLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, WernerWengleLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(LawOfTheWall, WernerWengleLawOfTheWall, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall
(
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    A_(constDict_.lookupOrAddDefault<scalar>("A", 8.3)),
    B_(constDict_.lookupOrAddDefault<scalar>("B", 0.14285714285714285))
{
    if (debug)
    {        
        printCoeffs();
    }
    
}

Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    WernerWengleLawOfTheWall(dict)
{
}


Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall
(
    const scalar A,
    const scalar B
)
:
    LawOfTheWall(),
    A_(A),
    B_(B)
{

    constDict_.add("A", A);
    constDict_.add("B", B);
    
    if (debug)
    {        
        printCoeffs();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::WernerWengleLawOfTheWall::printCoeffs() const
{
    Info<< nl << "WernerWengle law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "A" << indent << A_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::WernerWengleLawOfTheWall::value
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{  
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");

    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    scalar y = sampler.h()[index];
    return value(u, y, uTau, nu);
}


Foam::scalar Foam::WernerWengleLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{  
    scalar uPlus = u/uTau;
    scalar yPlus = y*uTau/nu;
    scalar yPlusM = pow(A_, 1/(1-B_));
    
    if (yPlus <= yPlusM)
    {
        return uPlus - yPlus;
    }
    else
    {
        return uPlus - A_*pow(yPlus, B_);
    }
}


Foam::scalar Foam::WernerWengleLawOfTheWall::derivative
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu        
) const
{
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");

    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    scalar y = sampler.h()[index];
    return derivative(u, y, uTau, nu);
}


Foam::scalar Foam::WernerWengleLawOfTheWall::derivative
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu        
) const
{
    scalar yPlus = y*uTau/nu;
    scalar yPlusM = pow(A_, 1/(1-B_));

    if (yPlus <= yPlusM)
    {
        return -u/sqr(uTau) - y/nu;
    }
    else
    {
        return -u/sqr(uTau) - A_*B_*pow(y/nu, B_)*pow(uTau, B_-1);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
