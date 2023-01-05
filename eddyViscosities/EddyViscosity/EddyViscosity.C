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

#include "EddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(EddyViscosity, 0);
    defineRunTimeSelectionTable(EddyViscosity, Dictionary);
    defineRunTimeSelectionTable(EddyViscosity, TypeAndDictionary);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  
    
Foam::autoPtr<Foam::EddyViscosity> Foam::EddyViscosity::New 
(
    const dictionary & dict
)
{
    word modelName(dict.lookup("type"));
    
    auto cstrIter =
    DictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "EddyViscosity::New(const dictionary&)"
        )   << "Unknown EddyViscosity type "
            << modelName << nl << nl
            << "Valid EddyViscosity types are :" << nl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    
    dictionary temp(dict);
    temp.remove("type");
   
    return cstrIter()(temp);  
}
 
Foam::autoPtr<Foam::EddyViscosity> Foam::EddyViscosity::New 
(
    const word & modelName,
    const dictionary & dict
)
{
    
    auto cstrIter =
    TypeAndDictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == TypeAndDictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "EddyViscosity::New(const word & modelName, const dictionary&)"
        )   << "Unknown EddyViscosity type "
            << modelName << nl << nl
            << "Valid EddyViscosity types are :" << nl
            << TypeAndDictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(modelName, dict);  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EddyViscosity::write(Ostream & os) const
{
    
    auto keys = constDict_.keys();
    label dictSize = constDict_.keys().size();
    
    os.writeKeyword("EddyViscosity")
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
