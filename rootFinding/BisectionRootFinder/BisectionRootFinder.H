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

    Controlled by the maximum number of iterations. The bisection method is
    meant to be used with equations where the Jacobian is not well defined.
    The root must be bracketed by the supplied lower and upper bounds.
    The convergence tolerance is fixed internally as a binary digit target for
    the root estimate.

    Usage
    \verbatim
    RootFinder
    {
        type       Bisection;
        maxIter    value; (default 30)
    }
    \endverbatim

Contributors/Copyright:
    2016-2026 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    BisectionRootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef BisectionRootFinder_H
#define BisectionRootFinder_H

#include "RootFinder.H"
#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BisectionRoot Declaration
\*---------------------------------------------------------------------------*/

class BisectionRootFinder : public RootFinder
{
private:

    //- Number of binary digits requested in the solution
    label getDigits_ =
        static_cast<int>(std::numeric_limits<scalar>::digits * 0.4);

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
