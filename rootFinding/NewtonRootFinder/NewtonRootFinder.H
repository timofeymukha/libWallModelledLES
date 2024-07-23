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
    Foam::NewtonRootFinder

@brief
    Root finder based on Newton's method.

    Controlled by two parameters: the max number of iterations and the relative
    error tolerance.

    Usage
    \verbatim
    RootFinder
    {
        type       Newton;
        maxIter    value; (default 30)
    }
    \endverbatim

Contributors/Copyright:
    2016-2019 Timofey Mukha
    2017      Saleh Rezaeiravesh

SourceFiles
    NewtonRootFinder.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonRootFinder_H
#define NewtonRootFinder_H

#include "RootFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NewtonRoot Declaration
\*---------------------------------------------------------------------------*/

class NewtonRootFinder : public RootFinder
{
private:
    //- Number of correct digits in the solution
    label getDigits_ =
        static_cast<int>(std::numeric_limits<scalar>::digits * 0.6);

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    TypeName("Newton");
#endif

    // Constructors

        //- Construct given name, function, its derivative, tolerance and
        //  maximum number of iterations
        NewtonRootFinder
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
        NewtonRootFinder
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
        NewtonRootFinder
        (
            const dictionary & dict
        )
        :
        RootFinder(dict)
        {}

        //- Copy constructor
        NewtonRootFinder(const NewtonRootFinder &) = default;

        //- Clone the object
        virtual autoPtr<RootFinder> clone() const
        {
            return autoPtr<RootFinder>
            (
                new NewtonRootFinder(*this)
            );
        }


    //- Destructor
        virtual ~NewtonRootFinder(){};

    // Member Functions

        //- Compute and return root
        std::pair<scalar, label> root
        (
            scalar guess,
            scalar lowerBound,
            scalar upperBound
        ) const;

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
