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
    NewtonRootFinder

Description
    Root finder based on Newton's method.
Author
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    NewtonRootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRootFinder_H
#define NewtonRootFinder_H

#include "RootFinder.H"
#include "scalar.H"
#include "word.H"
#include "typeInfo.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NewtonRoot Declaration
\*---------------------------------------------------------------------------*/

class NewtonRootFinder : public RootFinder
{

public:

    TypeName("Newton");
    
    // Constructors

        //- Construct given name, function, its derivative, tolerance and max number of iterations
        NewtonRootFinder
        (
            const word & rootFinderName,
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const scalar eps,
            const label maxIter
        )
        :
        RootFinder(rootFinderName, f, d, eps, maxIter)
        {}
        
        //- Construct given a function, its derivative, and dictionary
        NewtonRootFinder
        (
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary & dict
        )
        :
        RootFinder(f, d, dict)
        {}         
        
        
        //- Construct given  dictionary
        NewtonRootFinder
        (
            const dictionary & dict
        )
        :
        RootFinder(dict)
        {}         
       
        
    // Destructor - default
        virtual ~NewtonRootFinder(){};

    // Member Functions

        //- Return root
        scalar root(scalar x0) const;
        
        //- Write
        virtual void write(Ostream& os) const
        {
            RootFinder::write(os);
            os << decrIndent;
            os.writeKeyword("}") << endl;
        }

};
       



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
