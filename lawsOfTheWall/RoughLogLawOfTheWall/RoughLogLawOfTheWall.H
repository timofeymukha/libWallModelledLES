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
    Foam::RoughLogLawOfTheWall

@brief
    A log-law, in which the wall-normal distance is scaled by roughness
    length-scale.

    \f[
    u^+ =  1/\kappa \log \left( y/k_s \right) + B
    \f]

    Since the formula above provides an explicit value for the friction velocity
    Newton's method will converge in a single iteration.

    Usage:
    \verbatim
    Law
    {
        type      RoughLogLaw;
        kappa     value; (default 0.4)
        B         value; 
        ks        value;
    }
    \endverbatim

Contributors/Copyright:
    2019-2023 Timofey Mukha

SourceFiles
    RoughLogLawOfTheWall.C

\*---------------------------------------------------------------------------*/

#ifndef RoughLogLawOfTheWall_H
#define RoughLogLawOfTheWall_H

#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RoughLogLawOfTheWall Declaration
\*---------------------------------------------------------------------------*/

class RoughLogLawOfTheWall: public LawOfTheWall
{
    
    // Private data
        
        //- The kappa model constant
        scalar kappa_;
        
        //- The B model constant
        scalar B_;

        //- The roughness length-scale
        scalar ks_;
        
public:
    
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("RoughLogLaw");
#endif

    // Constructors
      
        //- Construct provided dictionary
        RoughLogLawOfTheWall(const dictionary &);
        
        //- Construct provided TypeName and dictionary
        RoughLogLawOfTheWall
        (
            const word & lawname,
            const dictionary &
        );

        //- Construct from model constants
        RoughLogLawOfTheWall
        (
            const scalar kappa,
            const scalar B,
            const scalar ks
        );

        //- Copy constructor
        RoughLogLawOfTheWall(const RoughLogLawOfTheWall &) = default;

        //- Assignment
        RoughLogLawOfTheWall & operator=(const RoughLogLawOfTheWall &) = default;

        //- Clone
        virtual autoPtr<LawOfTheWall> clone() const override
        {
            return autoPtr<LawOfTheWall>
            (
                new RoughLogLawOfTheWall(*this) 
            );
        }

    //- Destructor
        virtual ~RoughLogLawOfTheWall() {};

    // Member Functions
        
        scalar kappa() const
        {
            return kappa_;
        }
        
        scalar B() const
        {
            return B_;
        }

        scalar ks() const
        {
            return ks_;
        }
        
        //- Print the model coefficients to terminal
        virtual void printCoeffs() const override;

        //- Return the value of the implicit function defining the law
        virtual scalar value
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar value(scalar u, scalar y, scalar uTau, scalar nu) const;

                             
        //- Return the value of the  derivative of the implicit function
        //  defining the law.
        virtual scalar derivative
        (
            const SingleCellSampler & sampler,
            label index,
            scalar uTau,
            scalar nu
        ) const override;

        scalar derivative() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
