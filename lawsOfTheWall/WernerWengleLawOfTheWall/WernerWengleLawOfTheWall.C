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
    WernerWengleLawOfTheWall

Description
    The law of the wall proposed by WernerWengle.

Authors
    Timofey Mukha.  All rights reserved.

 \*---------------------------------------------------------------------------*/

#include "WernerWengleLawOfTheWall.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(WernerWengleLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, WernerWengleLawOfTheWall, Dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall()
:
    A_(8.3),
    B_(1./7)
{
    printCoeffs();
}


Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall
(
    scalar A,
    scalar B
)
:
    A_(A),
    B_(B)
{
    printCoeffs();
}

Foam::WernerWengleLawOfTheWall::WernerWengleLawOfTheWall
(
    const Foam::dictionary & dict
)
:
    LawOfTheWall(dict),
    A_(dict.lookupOrDefault<scalar>("A", 8.3)),
    B_(dict.lookupOrDefault<scalar>("B", 1./7))
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::WernerWengleLawOfTheWall::printCoeffs() const
{
    Info<< nl << "WernerWengle law of the wall" << nl;
    Info<< "{" << incrIndent << nl;
    Info<< indent << "A" << indent << A_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< "}" << nl << nl;
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
// ************************************************************************* //
}