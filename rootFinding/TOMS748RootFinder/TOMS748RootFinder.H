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
    Foam::TOMS748RootFinder

@brief
    Root finder wrapping Boost's toms748 algorithm.

    Controlled by on parameter: the max number of iterations.

    Usage
    \verbatim
    RootFinder
    {
        type       TOMS748;
        maxIter    value; (default 30)
    }
    \endverbatim

Contributors/Copyright:
    2023 Timofey Mukha

SourceFiles
    TOMS748RootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef TOMS748RootFinder_H
#define TOMS748RootFinder_H

#include "RootFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class TOMS748Root Declaration
\*---------------------------------------------------------------------------*/

class TOMS748RootFinder : public RootFinder
{
private:
    //- Number of correct digits in the solution
    label getDigits_ =
        static_cast<int>(std::numeric_limits<scalar>::digits - 5);

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("TOMS748");
#endif

    // Constructors

        //- Construct given name, function, its derivative, tolerance and
        //  maximum number of iterations
        TOMS748RootFinder
        (
            const word & rootFinderName,
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const label maxIter
        )
        :
        RootFinder(rootFinderName, f, d, maxIter)
        {}

        //- Construct given a function, its derivative, and dictionary
        TOMS748RootFinder
        (
            std::function<scalar(scalar)> f,
            std::function<scalar(scalar)> d,
            const dictionary & dict
        )
        :
        RootFinder(f, d, dict)
        {
        }


        //- Construct given  dictionary
        TOMS748RootFinder
        (
            const dictionary & dict
        )
        :
        RootFinder(dict)
        {}

        //- Copy constructor
        TOMS748RootFinder(const TOMS748RootFinder &) = default;

        //- Clone the object
        virtual autoPtr<RootFinder> clone() const
        {
            return autoPtr<RootFinder>
            (
                new TOMS748RootFinder(*this)
            );
        }


    //- Destructor
        virtual ~TOMS748RootFinder(){};

    // Member Functions

        //- Compute and return root
        scalar root(scalar guess, scalar lowerBound, scalar upperBound) const;

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
