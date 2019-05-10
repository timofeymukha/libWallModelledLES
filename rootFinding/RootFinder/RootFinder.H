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
    Foam::RootFinder

Description
    Base abstract class for finding roots of non-linear algebraic equations.

    In the context of wall modelling the latter are provided by laws of the
    wall. The root finders are therefore used in conjuction with the LOTW
    wall model
 
Contributors/Copyright:
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh


\*---------------------------------------------------------------------------*/

#ifndef RootFinder_H
#define RootFinder_H

#include "dictionary.H"
#include "refCount.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Dummy function to initialize the root finder
scalar rootFinderDummy(scalar);


/*---------------------------------------------------------------------------*\
                          Class RootFinder Declaration
\*---------------------------------------------------------------------------*/
   
class RootFinder : public refCount
{

protected:     

    //- Function defining the equation to be solved 
    std::function<scalar(scalar)> f_;
    
    //- Function defining the derivative of the equation to be solved
    std::function<scalar(scalar)> d_;
    
     //- Target relative error
    const scalar eps_;
        
    //- Maximum number of iterations
    const label maxIter_;


public:

    // Static data members
        TypeName ("RootFinder");

    // Constructors

        RootFinder 
        (        
            const word& rootFinderName,
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const scalar eps,
            const label maxIter
        )
        :
        f_(f),
        d_(d),
        eps_(eps),
        maxIter_(maxIter)
        {};
             
        RootFinder 
        (        
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary& dict
        )
        :
        f_(f),
        d_(d),
        eps_(dict.lookupOrDefault<scalar>("eps", 1e-3)),
        maxIter_(dict.lookupOrDefault<label>("maxIter", 30))
        {};
        
        RootFinder 
        (        
            const dictionary& dict
        )
        :
        f_(rootFinderDummy),
        d_(rootFinderDummy),
        eps_(dict.lookupOrDefault<scalar>("eps", 1e-3)),
        maxIter_(dict.lookupOrDefault<label>("maxIter", 30))
        {};

        RootFinder(const RootFinder &) = default;


    //- Destructor
        virtual ~RootFinder() {};
    
    // Selectors
        static autoPtr<RootFinder> New 
        (
            const word& rootFinderName,
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const scalar eps,
            const label maxIter
        );

        static autoPtr<RootFinder> New 
        (
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary & dict
        );

        static autoPtr<RootFinder> New 
        (
            const dictionary & dict
        );
   
    // Member Functions

        //- Return root
        virtual scalar root(scalar) const = 0;
        
        //- Set the implicit function defining the equation
        void setFunction(std::function<scalar(scalar)> f)
        {
            f_ = f;
        }
        
        //- Set the implicit function defining the derivative
        virtual void setDerivative(std::function<scalar(scalar)> d)
        {
            d_ = d;
        }
        
        //- Return eps_
        scalar eps() const
        {
            return eps_;
        }
                
        //- Return maxIter_
        label maxIter() const
        {
            return maxIter_;
        }
        
        //- Return the function f_
        std::function<scalar(scalar)> f() const
        {
            return f_;
        }
        
        //- Return the function d_
        std::function<scalar(scalar)> d() const
        {
            return d_;
        }

        //- Compute the value of f
        scalar f(scalar v) const
        {
            return f_(v);
        }

        //- Compute the value of d
        scalar d(scalar v) const
        {
            return d_(v);
        }
        
        virtual void write(Foam::Ostream& os) const;

    // RTS tables 
        
        // RTS table "Word"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            RootFinder,
            Word,
            (
                const word& rootFinderName,
                std::function<scalar(scalar)> f,
                std::function<scalar(scalar)> d,
                const scalar eps,
                const label maxIter
            ),
            (rootFinderName, f, d, eps, maxIter)
        );

        // RTS table "Dictionary"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            RootFinder,
            Dictionary,
            (
                std::function<scalar(scalar)> f,
                std::function<scalar(scalar)> d,
                const dictionary & dict
            ),
            (f, d, dict)
        );

        // RTS table "DictionaryOnly"    
        declareRunTimeSelectionTable
        (
            autoPtr,
            RootFinder,
            DictionaryOnly,
            (
                const dictionary & dict
            ),
            (dict)
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
