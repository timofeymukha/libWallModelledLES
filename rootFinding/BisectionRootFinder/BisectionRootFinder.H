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
    BisectionRootFinder

Description
    Root finder based on the bisection method.
Author
    Timofey Mukha. All rights reserved.

SourceFiles
    BisectionRootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef BisectionRootFinder_H
#define BisectionRootFinder_H

#include "RootFinder.H"
#include "scalar.H"
#include "word.H"
#include "typeInfo.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BisectionRoot Declaration
\*---------------------------------------------------------------------------*/

class BisectionRootFinder : public RootFinder
{
    // Private data
    
    //- The bracket defining the interval where the root will be searched
    //  [1/bracket*guess; bracket*guess]
    scalar bracket_;

public:

    TypeName("Bisection");
    
    // Constructors

        //- Construct given a function and its derivative
        BisectionRootFinder
        (
            const word & rootFinderName,
            std::function<scalar(scalar)>,
            std::function<scalar(scalar)>,
            const scalar eps,
            const label maxIter
        );

        //- Construct given a function, its derivative and dict
        BisectionRootFinder
        (
            std::function<scalar(scalar)>,
            std::function<scalar(scalar)>,
            const dictionary & dict
        );
        
        //- Construct given dicttionary
        BisectionRootFinder
        (
            const dictionary & dict
        );
        
    // Destructor - default
        virtual ~BisectionRootFinder(){};

    // Member Functions

        //- return the bracket
        scalar bracket() const
        {
            return bracket_;
        }
        
        //- Return root
        scalar root(scalar x0) const;
        
        
        virtual void write(Ostream& os) const
        {
            RootFinder::write(os);
            os.writeKeyword("bracket") << bracket_ << token::END_STATEMENT
                                       << endl;
            os << decrIndent;
            os.writeKeyword("}") << endl;
        }
        
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
