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
    EddyViscosity

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
#include "EddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EddyViscosity, 0);
    defineRunTimeSelectionTable(EddyViscosity, Dictionary);
    defineRunTimeSelectionTable(EddyViscosity, TypeAndDictionary);
    
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  
    
autoPtr<EddyViscosity> EddyViscosity::New 
(
    const dictionary & dict
)
{
//    word lawName(dict.lookup("type"));
    word modelName(dict.lookup("type"));
    
    DictionaryConstructorTable::iterator cstrIter =
//    DictionaryConstructorTablePtr_->find(lawName);
    DictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "EddyViscosity::New(const dictionary&)"
        ) << "Unknown EddyViscosity type "
//        << lawName << nl << nl
        << modelName << nl << nl
        << "Valid EddyViscosity types are :" << endl
        << DictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }
    
    dictionary temp(dict);
    temp.remove("type");
    
    Info << "Constructing EddyVisocity E1" <<nl;

    return cstrIter()(temp);  
}
 
autoPtr<EddyViscosity> EddyViscosity::New 
(
    const word & modelName,
    const dictionary & dict
)
{
    
    TypeAndDictionaryConstructorTable::iterator cstrIter =
    TypeAndDictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == TypeAndDictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "EddyViscosity::New(const word & modelName, const dictionary&)"
        ) << "Unknown EddyViscosity type "
        << modelName << nl << nl
        << "Valid EddyViscosity types are :" << endl
        << TypeAndDictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return cstrIter()(modelName, dict);  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void EddyViscosity::write(Ostream & os) const
{
    
    auto keys = constDict_.keys();
    label dictSize = constDict_.keys().size();
    
//    os.writeKeyword("Law") << endl;
    os.writeKeyword("EddyViscosity") << endl;
    os.writeKeyword("{") << incrIndent << endl;
    os.writeKeyword("type") << type() << token::END_STATEMENT << endl;
   
    for (int i=0; i<dictSize; i++)
    {
        os.writeKeyword(keys[i]) << constDict_[keys[i]][0] 
                                 << token::END_STATEMENT  << endl;
    }
   
    os << decrIndent;
    os.writeKeyword("}") << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
