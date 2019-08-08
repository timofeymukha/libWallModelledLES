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

#include "SampledField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SampledField, 0);
    defineRunTimeSelectionTable(SampledField, FvPatch);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::SampledField> Foam::SampledField::New 
(
    const word & fieldName,
    const fvPatch & p
)
{
    FvPatchConstructorTable::iterator cstrIter =
        FvPatchConstructorTablePtr_->find(fieldName);

    if (cstrIter == FvPatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "SampleDField::New(const word&, const fvPatch & p"
            
        )   << "Unknown SampledField type "
            << fieldName << nl << nl
            << "Valid SampledField types are :" << nl
            << FvPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(fieldName, p);
}  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //