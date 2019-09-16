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

#include "IntegratedWernerWengleLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(IntegratedWernerWengleLawOfTheWall, 0);
    addToRunTimeSelectionTable
    (
            LawOfTheWall,
            IntegratedWernerWengleLawOfTheWall,
            Dictionary);
    addToRunTimeSelectionTable
    (
            LawOfTheWall,
            IntegratedWernerWengleLawOfTheWall,
            TypeAndDictionary
    );
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IntegratedWernerWengleLawOfTheWall::IntegratedWernerWengleLawOfTheWall
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

Foam::IntegratedWernerWengleLawOfTheWall::IntegratedWernerWengleLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    IntegratedWernerWengleLawOfTheWall(dict)
{}



Foam::IntegratedWernerWengleLawOfTheWall::IntegratedWernerWengleLawOfTheWall
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

void Foam::IntegratedWernerWengleLawOfTheWall::printCoeffs() const
{
    Info<< nl << "IntegratedWernerWengle law of the wall" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "A" << indent << A_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}


Foam::scalar Foam::IntegratedWernerWengleLawOfTheWall::value
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");

    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    
    scalar h = sampler.h()[index]; 
    //scalar h1 = h - sampler.lengthList()[index]/2;
    scalar h2 = h + sampler.lengthList()[index]/2;
    return value(u, h2, uTau, nu);
}

Foam::scalar Foam::IntegratedWernerWengleLawOfTheWall::value
(
    scalar u,
    scalar h,
    scalar uTau,
    scalar nu
) const
{
    return uTau - pow((1+B_)/A_*pow(nu/h, B_)*u + 
                      (1-B_)/2*pow(A_, (1+B_)/(1-B_))*pow(nu/h, B_+1),
                      1.0/(1+B_));
}


Foam::scalar Foam::IntegratedWernerWengleLawOfTheWall::derivative
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu        
) const
{
    return derivative();
}

Foam::scalar Foam::IntegratedWernerWengleLawOfTheWall::derivative() const
{
    return 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
