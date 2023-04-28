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

#include "RoughLogLawOfTheWall.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "scalarListIOList.H"
#include "SingleCellSampler.H"
#include "codeRules.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(RoughLogLawOfTheWall, 0);
    addToRunTimeSelectionTable(LawOfTheWall, RoughLogLawOfTheWall, Dictionary);
    addToRunTimeSelectionTable(LawOfTheWall, RoughLogLawOfTheWall, TypeAndDictionary);

}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RoughLogLawOfTheWall::RoughLogLawOfTheWall
(
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
#ifdef FOAM_DICTIONARY_NO_GET
#ifdef FOAM_DICTIONARY_HAS_LOOKUP
    B_(constDict_.lookup<scalar>("B")),
    ks_(constDict_.lookup<scalar>("ks"))
#else
    B_(constDict_.lookupType<scalar>("B")),
    ks_(constDict_.lookupType<scalar>("ks"))
#endif
#else
    B_(constDict_.get<scalar>("B")),
    ks_(constDict_.get<scalar>("ks"))
#endif
{
    if (debug)
    {        
        printCoeffs();
    }
}


Foam::RoughLogLawOfTheWall::RoughLogLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    RoughLogLawOfTheWall(dict)
{
}


Foam::RoughLogLawOfTheWall::RoughLogLawOfTheWall
(
    const scalar kappa,
    const scalar B,
    const scalar ks
)
:
    LawOfTheWall(),
    kappa_(kappa),
    B_(B),
    ks_(ks)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B", B);
    constDict_.add("ks", ks_);

    if (debug)
    {        
        printCoeffs();
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RoughLogLawOfTheWall::printCoeffs() const
{
    Info<< nl << "RoughLogLaw law of the wall" << nl;     
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "kappa" << indent << kappa_ << nl;
    Info<< indent << "B" << indent <<  B_ << nl;
    Info<< indent << "ks" << indent <<  ks_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}

Foam::scalar Foam::RoughLogLawOfTheWall::value
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
    return  value(u, y, uTau, nu);
}

Foam::scalar Foam::RoughLogLawOfTheWall::value
(
    scalar u,
    scalar y,
    scalar uTau,
    scalar nu
) const
{
    scalar denom = 1/kappa_*log(y/ks_) + B_;
    return uTau - u/denom;
}

Foam::scalar Foam::RoughLogLawOfTheWall::derivative
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu        
) const
{
    return  1;
}

Foam::scalar Foam::RoughLogLawOfTheWall::derivative() const
{
    return 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

