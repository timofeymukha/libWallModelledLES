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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SampledField::projectVectors(scalarListList & field) const
{
    const tmp<vectorField> tfaceNormals = patch_.nf();
    const vectorField faceNormals = tfaceNormals();

    
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

void Foam::SampledField::projectVectors(scalarListListList & field) const
{
    const tmp<vectorField> tfaceNormals = patch_.nf();
    const vectorField faceNormals = tfaceNormals();

    
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

void Foam::SampledField::projectVectors(vectorField & field) const
{
    const tmp<vectorField> tfaceNormals = patch_.nf();
    const vectorField faceNormals = tfaceNormals();

    
    forAll(field, i)
    {   
        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(field[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        field[i] -= normal;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
