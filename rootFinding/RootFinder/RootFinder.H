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
    RootFinder

Description
    Base class for finding roots of non-linear algebraic equations.
Author
    Timofey Mukha. All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef RootFinder_H
#define RootFinder_H

#include "scalar.H"
#include "dictionary.H"
#include "word.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RootFinder Declaration
\*---------------------------------------------------------------------------*/

class RootFinder
{

protected:     

    //- Pointer to the function
    std::function<scalar(scalar)> f_;
    
     //- Tolerance
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
        eps_(dict.lookupOrDefault<scalar>("eps", 1e-2)),
        maxIter_(dict.lookupOrDefault<label>("maxIter", 15))
        {};

    // Destructor
        virtual ~RootFinder() {};

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
    )
    
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
    )
        
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
   
    // Member Functions

        //- Return root
        virtual scalar root(scalar guess) const = 0;
        
        //- Set the implicit function defining the equation
        void setFunction(std::function<scalar(scalar)> f)
        {
            f_ = f;
        }
        
        //- Set the implicit function defining the derivative
        virtual void setDerivative(std::function<scalar(scalar)> d){}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
