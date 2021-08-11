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

#include "DupratEddyViscosity.H"
#include "dictionary.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"
#include "SampledPGradField.H"
#include "scalarListIOList.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(DupratEddyViscosity, 0);
    addToRunTimeSelectionTable(EddyViscosity, DupratEddyViscosity, Dictionary);
    addToRunTimeSelectionTable(EddyViscosity, DupratEddyViscosity, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DupratEddyViscosity::DupratEddyViscosity
(
    const dictionary & dict
)
:
    EddyViscosity(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    APlus_(constDict_.lookupOrAddDefault<scalar>("APlus", 18)),
    beta_(constDict_.lookupOrAddDefault<scalar>("beta", 0.78))
{
    if (debug)
    {        
        printCoeffs();
    }
    
}

Foam::DupratEddyViscosity::DupratEddyViscosity
(
    const word & modelName,
    const dictionary & dict
)
:
    DupratEddyViscosity(dict)
{
}

Foam::DupratEddyViscosity::DupratEddyViscosity
(
    const scalar kappa,
    const scalar APlus,
    const scalar beta
)
:
    EddyViscosity(),
    kappa_(kappa),
    APlus_(APlus),
    beta_(beta)
{
    constDict_.add("kappa", kappa);
    constDict_.add("APlus", APlus);
    constDict_.add("beta", beta);

    if (debug)
    {        
        printCoeffs();
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DupratEddyViscosity::addFieldsToSampler(Sampler & sampler) const
{
    if (debug)
    {
        Info<< "Duprat eddy viscosity: Adding pressure gradient to sampler"
            << nl;
    }
    sampler.addField(new SampledPGradField(sampler.patch()));
}

void Foam::DupratEddyViscosity::printCoeffs() const
{
    Info<< nl << "Duprat eddy viscosity model" << nl;
    Info<< token::BEGIN_BLOCK << incrIndent << nl;
    Info<< indent << "APlus" << indent << APlus_ << nl;
    Info<< indent << "kappa" << indent <<  kappa_ << nl;
    Info<< indent << "beta" << indent <<  beta_ << nl;
    Info<< token::END_BLOCK << nl << nl;
}

Foam::scalarList Foam::DupratEddyViscosity::value
(
    const SingleCellSampler & sampler,
    const label index,
    const scalarList & y,
    const scalar uTau,
    const scalar nu
) const
{  
    const scalarListList & pGrad =
        sampler.db().lookupObject<scalarListIOList>("pGrad");

    scalarList pGradI = pGrad[index];
    scalar magPGrad = mag(vector(pGradI[0], pGradI[1], pGradI[2]));

    return value(y, magPGrad, uTau, nu);
}

Foam::scalarList Foam::DupratEddyViscosity::value
(
    const scalarList & y,
    const scalar magPGrad,
    const scalar uTau,
    const scalar nu
) const
{  

    const scalar uP = pow(nu*magPGrad, 1./3);
    const scalar uTauP = sqrt(sqr(uTau) + sqr(uP));
    const scalar alpha = sqr(uTau)/sqr(uTauP);
        
    const scalarList yStar = y*uTauP/nu;
    
    scalarList values(y.size(), 0.0);
     
    forAll(values, i)
    {
        values[i] = nu*kappa_*yStar[i]*
                    pow(alpha + yStar[i]*pow(1 - alpha, 1.5), beta_)*
                    sqr(1 - exp(-yStar[i]/(1 + APlus_*pow(alpha, 3))));
    }
    return values;
}

// ************************************************************************* //

