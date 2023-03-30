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

\*---------------------------------------------------------------------------*/

#include "RootFinder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(RootFinder, 0);
    defineRunTimeSelectionTable(RootFinder, Word);
    defineRunTimeSelectionTable(RootFinder, Dictionary);
    defineRunTimeSelectionTable(RootFinder, DictionaryOnly);
}
#endif


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::RootFinder> Foam::RootFinder::New 
(
    const word & rootFinderName,
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const scalar eps,
    const label maxIter
)
{
    auto cstrIter =
    WordConstructorTablePtr_->find(rootFinderName);

    if (cstrIter == WordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RootFinder::New(const word&, std::function<scalar(scalar)> f,"
                "std::function<scalar(scalar)> d, const scalar eps,"
                "const label maxIter)"
        )   << "Unknown RootFinder type "
            << rootFinderName << nl << nl
            << "Valid RootFinder types are :" << nl
            << WordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(rootFinderName, f, d, eps, maxIter);
}  


Foam::autoPtr<Foam::RootFinder> Foam::RootFinder::New 
(
    std::function<scalar(scalar)> f,
    std::function<scalar(scalar)> d,
    const Foam::dictionary & dict
)
{
    Foam::word rootFinderName(dict.lookup("type"));
    //const entry* ePtr = dict.lookupEntryPtr("type", false, false);
    //word rootFinderName(ePtr->stream());
    //ePtr->stream() >> rootFinderName;

    auto cstrIter =
    DictionaryConstructorTablePtr_->find(rootFinderName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RootFinder::New(std::function<scalar(scalar)> f,"
                "std::function<scalar(scalar)> d, const dictionary & dict)"
        )   << "Unknown RootFinder type "
            << rootFinderName << nl << nl
            << "Valid RootFinder types are :" << nl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(f, d, dict);
}


Foam::autoPtr<Foam::RootFinder> Foam::RootFinder::New 
(
    const Foam::dictionary & dict
)
{
    Foam::word rootFinderName(dict.lookup("type"));

    auto cstrIter =
        DictionaryOnlyConstructorTablePtr_->find(rootFinderName);

    if (cstrIter == DictionaryOnlyConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RootFinder::New(const dictionary & dict)"
        )   << "Unknown RootFinder type "
            << rootFinderName << nl << nl
            << "Valid RootFinder types are :" << nl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RootFinder::write(Foam::Ostream & os) const
{
    os.writeKeyword("RootFinder") << nl;
    os.writeKeyword("{")
        << incrIndent << nl;
    os.writeKeyword("type")
        << type() << token::END_STATEMENT << nl;
    os.writeKeyword("eps")
        << eps_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxIter")
        << maxIter_ << token::END_STATEMENT << nl; 
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
