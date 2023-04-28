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
    Foam::BisectionRootFinder

@brief
    Root finder based on the bisection method. 
 
   Controlled by three paramters:
    the max number of iterations, the relative error tolerance and the bracket.
    The latter is used to define the initial size of the search interval in
    conjuction with an initial guess of the root value. If x0 is the guess then
    the search interval is from x0/bracket to bracket*x0.

    The bisection method is meant to be used with equations where the Jacobian
    is not well defined. Currently, however, it appears to be problematic to
    get a good initial interval size in the beginning of the simulation.

    Usage
    \verbatim
    RootFinder
    {
        type       Bisection;
        maxIter    value; (default 10)
        eps        value; (default 0.01)
        bracket    value; (default 3)
    }
    \endverbatim

Contributors/Copyright:
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    BisectionRootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef BisectionRootFinder_H
#define BisectionRootFinder_H

#include "RootFinder.H"

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

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Bisection");
#endif
    
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

        //- Construct given a function, its derivative and dictionary
        BisectionRootFinder
        (
            std::function<scalar(scalar)>,
            std::function<scalar(scalar)>,
            const dictionary & dict
        );
        
        //- Construct given dictionary
        BisectionRootFinder
        (
            const dictionary & dict
        );

        BisectionRootFinder(const BisectionRootFinder &) = default;

        //- Clone the object
        virtual autoPtr<RootFinder> clone() const
        {
            return autoPtr<RootFinder>
            (
                new BisectionRootFinder(*this)
            );
        }

    // Destructor - default
        virtual ~BisectionRootFinder(){};

    // Member Functions

        //- return the bracket
        scalar bracket() const
        {
            return bracket_;
        }
        
        //- Return root
        scalar root(scalar guess) const;
        
        //- Write parameters to stream
        virtual void write(Ostream& os) const
        {
            RootFinder::write(os);
            os.writeKeyword("bracket") 
                << bracket_ << token::END_STATEMENT << nl;   
            os  << decrIndent;
            os.writeKeyword("}")
                << endl;
        }       
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
