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
#include "scalarListIOList.H"
#include "SingleCellSampler.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
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
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IntegratedReichardtLawOfTheWall::IntegratedReichardtLawOfTheWall
(
    const dictionary & dict
)
:
    LawOfTheWall(dict),
    kappa_(constDict_.lookupOrAddDefault<scalar>("kappa", 0.4)),
    B1_(constDict_.lookupOrAddDefault<scalar>("B1", 11)),
    B2_(constDict_.lookupOrAddDefault<scalar>("B2", 3)),
    C_(constDict_.lookupOrAddDefault<scalar>("C", 7.8))
{
    if (debug)
    {        
        printCoeffs();
    }
    
}

Foam::IntegratedReichardtLawOfTheWall::IntegratedReichardtLawOfTheWall
(
    const word & lawName,
    const dictionary & dict
)
:
    IntegratedReichardtLawOfTheWall(dict)
{
}


Foam::IntegratedReichardtLawOfTheWall::IntegratedReichardtLawOfTheWall
(
    const scalar kappa,
    const scalar B1,
    const scalar B2,
    const scalar C
)
:
    LawOfTheWall(),
    kappa_(kappa),
    B1_(B1),
    B2_(B2),
    C_(C)
{
    constDict_.add("kappa", kappa);
    constDict_.add("B1", B1);
    constDict_.add("B2", B2);
    constDict_.add("C", C);

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
 
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{  
    const scalarListIOList & U = sampler.db().lookupObject<scalarListIOList>("U");
    scalar u = mag(vector(U[index][0], U[index][1], U[index][2]));
    
    scalar h = sampler.h()[index];
    // !!!!!
    scalar h1 = mag(h - sampler.lengthList()[index]/2);
    scalar h2 = h + sampler.lengthList()[index]/2;

    return value(u, h1, h2, uTau, nu);
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::valueMulticell
(
 
    const MultiCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu
) const
{  
    const scalarListList & U =
        sampler.db().lookupObject<scalarListListIOList>("U")[index];
    
    
    scalarList h = sampler.h()[index];
    scalarList l = sampler.lengthList()[index];

    // Compute cell-length weighted mean of u across sampling cells
    scalar uMean = 0;
    
    for(int i=0; i < h.size(); i++)
    {
        uMean += l[i]*mag(vector(U[i][0], U[i][1], U[i][2]));
    }

    scalar h1 = mag(h[0] - l[0]/2);
    scalar h2 = h[h.size()-1] + l[h.size()-1]/2;
    
    uMean = uMean/(h2 - h1);

    return value(uMean, h1, h2, uTau, nu);
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::value
(
    scalar u,
    scalar h1,
    scalar h2,
    scalar uTau,
    scalar nu
) const
{   

    
    return u*(h2 - h1) - (logTerm(h2, uTau, nu) - logTerm(h1, uTau, nu) +
                          expTerm(h2, uTau, nu) - expTerm(h1, uTau, nu));
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::derivative
(
    const SingleCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu        
) const
{ 
    scalar h = sampler.h()[index];
    
    scalar h1 = mag(h - sampler.lengthList()[index]/2);
    scalar h2 = h + sampler.lengthList()[index]/2;
        
    return derivative(h1, h2, uTau, nu);
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::derivativeMulticell
(
    const MultiCellSampler & sampler,
    label index,
    scalar uTau,
    scalar nu        
) const
{ 
    scalarList h = sampler.h()[index];
    scalar h1 = mag(h[0] - sampler.lengthList()[index][0]/2);
    scalar h2 = h[h.size()-1] +
                sampler.lengthList()[index][h.size()-1]/2;

    return derivative(h1, h2, uTau, nu);
}


Foam::scalar Foam::IntegratedReichardtLawOfTheWall::derivative
(
    scalar h1,
    scalar h2,
    scalar uTau,
    scalar nu        
) const
{ 
    return -(logTermDerivative(h2, uTau, nu) - logTermDerivative(h1, uTau, nu) +
             expTermDerivative(h2, uTau, nu) - expTermDerivative(h1, uTau, nu));
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
