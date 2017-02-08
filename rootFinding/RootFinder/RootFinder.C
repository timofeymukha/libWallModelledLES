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


#include "scalar.H"
#include "word.H"
#include "dictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "RootFinder.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RootFinder, 0);
    defineRunTimeSelectionTable(RootFinder, Word);
    defineRunTimeSelectionTable(RootFinder, Dictionary);


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

autoPtr<RootFinder> RootFinder::New 
(
    const word & rootFinderName,
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const scalar eps,
    const label maxIter
)
{
    WordConstructorTable::iterator cstrIter =
    WordConstructorTablePtr_->find(rootFinderName);

    if (cstrIter == WordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RootFinder::New(const word&, std::function<scalar(scalar)> f,"
                "std::function<scalar(scalar)> d, const scalar eps"
                "const label maxIter)"
        ) << "Unknown RootFinder type "
        << rootFinderName << nl << nl
        << "Valid RootFinder types are :" << endl
        << WordConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return cstrIter()(rootFinderName, f, d, eps, maxIter);
}  

autoPtr<RootFinder> RootFinder::New 
(
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const dictionary & dict
)
{
    word rootFinderName = dict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
    DictionaryConstructorTablePtr_->find(rootFinderName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RootFinder::New(std::function<scalar(scalar)> f,"
                "std::function<scalar(scalar)> d, const dictionary & dict)"
        ) << "Unknown RootFinder type "
        << rootFinderName << nl << nl
        << "Valid RootFinder types are :" << endl
        << DictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return cstrIter()(f, d, dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam