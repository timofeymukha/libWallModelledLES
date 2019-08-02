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

#include "VanDriestEddyViscosity.H"
#include "addToRunTimeSelectionTable.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(VanDriestEddyViscosity, 0);
    addToRunTimeSelectionTable(EddyViscosity, VanDriestEddyViscosity, Dictionary);
    addToRunTimeSelectionTable(EddyViscosity, VanDriestEddyViscosity, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VanDriestEddyViscosity::VanDriestEddyViscosity
(
    const dictionary & dict
)
:
    EddyViscosity(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    APlus_(constDict_.lookupOrAddDefault<scalar>("APlus", 18))
{

    if (debug)
    {        
        printCoeffs();
    }
    
}

Foam::VanDriestEddyViscosity::VanDriestEddyViscosity
(
    const word & modelName,
    const dictionary & dict
)
:
    VanDriestEddyViscosity(dict)
{
}

Foam::VanDriestEddyViscosity::VanDriestEddyViscosity
(
    const scalar kappa,
    const scalar APlus
)
:
    EddyViscosity(),
    kappa_(kappa),
    APlus_(APlus)
{
    constDict_.add("kappa", kappa);
    constDict_.add("APlus", APlus);

    if (debug)
    {        
        printCoeffs();
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VanDriestEddyViscosity::printCoeffs() const
{
    Info<< nl << "VanDriest eddy viscosity model" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "APlus" << indent << APlus_ << nl;
    Info<< indent << "kappa" << indent <<  kappa_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}

Foam::scalarList Foam::VanDriestEddyViscosity::value
(
    const SingleCellSampler & sampler,
    const label index, 
    const scalarList & y,
    const scalar uTau,
    const scalar nu
) const
{  
    return value(y, uTau, nu);
}


Foam::scalarList Foam::VanDriestEddyViscosity::value
(
    const scalarList & y,
    const scalar uTau,
    const scalar nu
) const
{  
    const scalarList yPlus = y*uTau/nu;
    
    scalarList values(y.size(), 0.0);
     
    forAll(values, i)
    {
        values[i] = kappa_*uTau*y[i]*sqr(1 - exp(-yPlus[i]/APlus_));
    }
    return values;
}

// ************************************************************************* //

