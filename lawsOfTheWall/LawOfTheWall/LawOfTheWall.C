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
    LawOfTheWall

Description
    Base class for laws of the wall.

Authors
    Timofey Mukha.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "dictionary.H"
#include "word.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "LawOfTheWall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LawOfTheWall, 0);
    defineRunTimeSelectionTable(LawOfTheWall, Dictionary);
    defineRunTimeSelectionTable(LawOfTheWall, TypeAndDictionary);
    
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  
    
autoPtr<LawOfTheWall> LawOfTheWall::New 
(
    const dictionary & dict
)
{
    word lawName(dict.lookup("type"));
    
    DictionaryConstructorTable::iterator cstrIter =
    DictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LawOfTheWall::New(const dictionary&)"
        ) << "Unknown LawOfTheWall type "
        << lawName << nl << nl
        << "Valid LawOfTheWall types are :" << endl
        << DictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return cstrIter()(dict);  
}
 
autoPtr<LawOfTheWall> LawOfTheWall::New 
(
    const word & lawName,
    const dictionary & dict
)
{
    
    TypeAndDictionaryConstructorTable::iterator cstrIter =
    TypeAndDictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == TypeAndDictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "LawOfTheWall::New(const word & lawName, const dictionary&)"
        ) << "Unknown LawOfTheWall type "
        << lawName << nl << nl
        << "Valid LawOfTheWall types are :" << endl
        << TypeAndDictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return cstrIter()(lawName, dict);  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LawOfTheWall::write(Ostream & os) const
{
    
    auto keys = constDict_.keys();
    label dictSize = constDict_.keys().size();
    
    os.writeKeyword("Law") << endl;
    os.writeKeyword("{") << endl;

    for (int i=0; i<dictSize; i++)
    {
        os.writeKeyword(keys[i]) <<  constDict_[keys[i]][0]<< token::END_STATEMENT  << endl;
    }
    if (!constDict_.found("type"))
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << endl;
        
    }
    os.writeKeyword("}") << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
