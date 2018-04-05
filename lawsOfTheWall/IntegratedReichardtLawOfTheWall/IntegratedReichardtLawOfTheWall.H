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

Class
    Foam::IntegratedReichardtLawOfTheWall

Description
    Model based on integrating the law of the wall proposed by Reichardt.

    \f[
    u^+ = \frac{1}{\kappa } \ln (1 + \kappa y^+) + C \left( 1 - \exp (-y^+/B_1)-
          \frac{y^+}{B_1} \exp(-y^+/B_2) \right)
    \f]   

    Usage:
    \verbatim
    Law
    {
        type      IntegratedReichardt;
        kappa     value; (default 0.4)
        B1        value; (default 11)
        B2        value; (default 3)
        C         value; (default 7.8)
    }
    \endverbatim

    Reference:
    \verbatim
        Reichardt, H. (1951).
        Vollstandige Darstellung der turbulenten Geschwindigkeitsverteilung in
        glatten Leitungen.
        Zeitschrift fur Angewandte Mathematik und Mechanik 31(7) (pp. 208-219).
    \endverbatim

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    IntegratedReichardtLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef IntegratedReichardtLawOfTheWall_H
#define IntegratedReichardtLawOfTheWall_H

#include "scalar.H"
#include "typeInfo.H"
#include "dictionary.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class WernerWengleLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class IntegratedReichardtLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa_;
        
        //- The B1 model constant
        scalar B1_;
        
        //- The B1 model constant
        scalar B2_;
        
        //- The C model constant
        scalar C_;
        
public:
    
    TypeName("IntegratedReichardt");

    // Constructors

        //- Construct provided dictionary and cell index list
        IntegratedReichardtLawOfTheWall
        (
            const dictionary &,
            const Sampler &
        );
        
        //- Construct provided TypeName and dictionary and cell index list
        IntegratedReichardtLawOfTheWall
        (
            const word & lawname,
            const dictionary &,
            const Sampler &    
        );

    // Destructor
        virtual ~IntegratedReichardtLawOfTheWall() {};

    // Member Functions
        
        //- Return the kappa constant
        scalar kappa() const
        {
            return kappa_;
        }
        
        //- Return the B1 constant
        scalar B1() const
        {
            return B1_;
        }
        
        //- Return the B2 constant
        scalar B2() const
        {
            return B2_;
        }
        
        //- Return the C constant
        scalar C() const
        {
            return C_;
        }
        
        //- Print info to terminal
        virtual void printCoeffs() const;

        //- Return the value of the implicit function defining the law
        virtual scalar value(scalar index, scalar uTau, scalar nu) const;
        
        //- Return the value of the derivative of the implicit function
        //  defining the law
        virtual scalar derivative(scalar index, scalar uTau, scalar nu) const;
        
        //- The log-term in the integrated law
        scalar logTerm(scalar y, scalar uTau, scalar nu) const;
        
        //- The exp-term in the integrated law
        scalar expTerm(scalar y, scalar uTau, scalar nu) const;
        
        //- The derivative of the log-term in the integrated law wrt uTau
        scalar logTermDerivative(scalar y, scalar uTau, scalar nu) const;
        
        //- The derivative of the exp-term in the integrated law wrt uTau
        scalar expTermDerivative(scalar y, scalar uTau, scalar nu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif