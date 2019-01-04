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

#include "SampledVelocityField.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::SampledVelocityField::sample
(
    Foam::scalarListList & sampledValues,
    const Foam::labelList & indexList
) const
{
    Info<< "Sampling velocity for patch " << patch_.name() << nl;
    
    const volVectorField & UField = db().lookupObject<volVectorField>("U");
    const vectorField & Uwall = UField.boundaryField()[patch().index()];
    
    vectorField sampledU(indexList.size());
    
    for (int i=0; i<indexList.size(); i++)
    {
        sampledU[i] = UField[indexList[i]] - Uwall[i]; 
        sampledValues[i] = List<scalar>(3);
        
        for (int j=0; j<3; j++)
        {
            sampledValues[i][j] = sampledU[i][j]; 
        }
    }
    projectVectors(sampledValues);
}


void Foam::SampledVelocityField::registerFields() const
{
    const objectRegistry & registry =
        db().subRegistry("wallModelSampling").subRegistry(patch_.name());

    scalarListList sampledU(patch().size());
        
    if (db().thisDb().foundObject<volVectorField>("U"))
    {
        forAll(sampledU, i)
        {
            sampledU[i] = scalarList(3, 0.0);
            //forAll(sampledU[i], j)
            //{
                //sampledU[i][j] = 0;
            //}
        }

        //projectVectors(sampledU);
    }


    db().thisDb().store
    (        
        new IOList<scalarList>
        (
            IOobject
            (
                "U",
                db().time().timeName(),
                registry, 
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sampledU
        )
    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
