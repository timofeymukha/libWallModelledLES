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
    Foam::PatchProjector

@brief
    Class for projecting fields and lists onto a patch.

Contributors/Copyright:
    2019 Timofey Mukha


SourceFiles
    helpers.C

\*---------------------------------------------------------------------------*/

#ifndef Helpers_H
#define Helpers_H

#include "fixedValueFvPatchFields.H"
#include "scalarListIOList.H"
#include "scalarListListIOList.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Helpers
{

    void projectOnPatch(tmp<vectorField> normals, vectorField & field);

    void projectOnPatch(tmp<vectorField> normals, scalarListListList & field);

    void projectOnPatch(tmp<vectorField> normals, scalarListList & field);


    tmp<scalarField> mag(const scalarListList & list);

    //- Convert listList to a field
    template<class Type>
    tmp<Field<Type> > listListToField(const scalarListList & list)
    {
        tmp<Field<Type> > tField
        (
            new Field<Type>(list.size(), pTraits<Type>::zero)
        );

#ifdef FOAM_NEW_TMP_RULES
        Field<Type> & field = tField.ref();
#else
        Field<Type> & field = tField();
#endif
        
        forAll(list, i)
        {
            Type element;
            forAll(list[i], j)
            {
                element[j] = list[i][j];
            }
            field[i] = element;
        }
        return tField;
    }

} // End namespace Foam

} // End namespace Helpers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

