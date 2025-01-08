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

   Controlled by two parameters:
    the max number of iterations and the relative error tolerance.
    The bisection method is meant to be used with equations where the Jacobian
    is not well defined.

    Usage
    \verbatim
    RootFinder
    {
        type       Bisection;
        maxIter    value; (default 10)
    }
    \endverbatim

Contributors/Copyright:
    2016-2025 Timofey Mukha
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

        //- Compute and return root
        std::pair<scalar, label> root
        (
            scalar guess,
            scalar lowerBound,
            scalar upperBound
        ) const;

        //- Write parameters to stream
        virtual void write(Ostream& os) const
        {
            RootFinder::write(os);
            os  << decrIndent;
            os.writeKeyword("}")
                << endl;
        }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
