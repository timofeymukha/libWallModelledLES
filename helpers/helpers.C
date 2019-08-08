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

#include "helpers.H"
#include "volFields.H"
#include "codeRules.H"


void Foam::Helpers::projectOnPatch
(
    tmp<vectorField> normals,
    vectorField & field
)
{
    const vectorField faceNormals = normals();

    forAll(field, i)
    {    
        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(field[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        field[i] -= normal;
    }
}

void Foam::Helpers::projectOnPatch
(
    tmp<vectorField> normals,
    scalarListList & field
)
{
    const vectorField faceNormals = normals();

    forAll(field, i)
    {
        vector fieldI(field[i][0], field[i][1], field[i][2]);

        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(fieldI & -faceNormals[i]);

        // Subtract normal component to get the parallel one
        fieldI -= normal;

        // Assign back to list
        List<scalar> fieldIList(3);
        fieldIList[0] = fieldI[0];
        fieldIList[1] = fieldI[1];
        fieldIList[2] = fieldI[2];
        field[i] = fieldIList;
    }
}

void Foam::Helpers::projectOnPatch
(
    tmp<vectorField> normals,
    scalarListListList & field
)
{
    const vectorField faceNormals = normals();

    forAll(field, i)
    {
        forAll(field[i], j)
        {
            vector fieldI(field[i][j][0], field[i][j][1], field[i][j][2]);

            // Normal component as dot product with (inwards) face normal
            vector normal = -faceNormals[i]*(fieldI & -faceNormals[i]);

            // Subtract normal component to get the parallel one
            fieldI -= normal;

            // Assign back to list
            List<scalar> fieldIList(3);
            fieldIList[0] = fieldI[0];
            fieldIList[1] = fieldI[1];
            fieldIList[2] = fieldI[2];
            field[i][j] = fieldIList;
        }
    }
}

Foam::tmp<Foam::scalarField> Foam::Helpers::mag(const scalarListList & list)
{
    tmp<scalarField> tField(new scalarField(list.size(), 0.0));

#ifdef FOAM_NEW_TMP_RULES
    scalarField & field = tField.ref();
#else
    scalarField & field = tField();
#endif
    
    forAll(list, i)
    {
        scalar element = 0;
        forAll(list[i], j)
        {
            element += sqr(list[i][j]);
        }
        field[i] = sqrt(element);
    }
    return tField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
