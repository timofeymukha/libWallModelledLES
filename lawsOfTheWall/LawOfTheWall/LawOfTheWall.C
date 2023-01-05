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

#include "LawOfTheWall.H"
#include "SingleCellSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(LawOfTheWall, 0);
    defineRunTimeSelectionTable(LawOfTheWall, Dictionary);
    defineRunTimeSelectionTable(LawOfTheWall, TypeAndDictionary);
}    
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  
    
Foam::autoPtr<Foam::LawOfTheWall> Foam::LawOfTheWall::New 
(
    const Foam::dictionary & dict
)
{
    Foam::word lawName(dict.lookup("type"));
    

    auto cstrIter =
        DictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LawOfTheWall::New(const dictionary&)"
        )   << "Unknown LawOfTheWall type "
            << lawName << nl << nl
            << "Valid LawOfTheWall types are :" << nl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    
    Foam::dictionary temp(dict);
    temp.remove("type");

    return cstrIter()(temp);  
}
 

Foam::autoPtr<Foam::LawOfTheWall> Foam::LawOfTheWall::New 
(
    const Foam::word & lawName,
    const Foam::dictionary & dict
)
{
    
    auto cstrIter =
        TypeAndDictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == TypeAndDictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LawOfTheWall::New(const word & lawName, const dictionary&)"
        )   << "Unknown LawOfTheWall type "
            << lawName << nl << nl
            << "Valid LawOfTheWall types are :" << nl
            << TypeAndDictionaryConstructorTablePtr_->sortedToc()
            << exit(Foam::FatalError);
    }

    return cstrIter()(lawName, dict);  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LawOfTheWall::write(Foam::Ostream & os) const
{
    
    auto keys = constDict_.keys();
    Foam::label dictSize = constDict_.keys().size();
    
    os.writeKeyword("Law") 
        << nl;
    os.writeKeyword("{")
        << incrIndent << nl;
    os.writeKeyword("type")
        << type() << token::END_STATEMENT << nl;
   
    for (int i=0; i<dictSize; i++)
    {
        os.writeKeyword(keys[i])
            << constDict_[keys[i]][0] << token::END_STATEMENT << nl;
    }
   
    os  << decrIndent;
    os.writeKeyword("}")
        << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
