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
#include "helpers.H"
#include "interpolation.H"
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SampledVelocityField, 0);
}
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SampledVelocityField::sample
(
    Foam::scalarListList & sampledValues,
    const Foam::labelList & indexList,
    const Foam::scalarField & h
) const
{
    if (debug)
    {
        Info<< "Sampling velocity for patch " << patch_.name() << nl;
    }

    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();

    const volVectorField & UField = mesh().lookupObject<volVectorField>("U");
    const vectorField & Uwall = UField.boundaryField()[patch().index()];

    vectorField sampledU(indexList.size());

    autoPtr<interpolation<vector> > interpolator;
    interpolator.operator=
    (
        interpolation<vector>::New(interpolationType(), UField)
    );

    for (int i=0; i<indexList.size(); i++)
    {

        point p = faceCentres[i] - h[i]*faceNormals[i];
        const vector interp = interpolator->interpolate(p, indexList[i]);
        sampledU[i] = interp - Uwall[i];
        scalarList temp(3, 0.0);
        
        for (int j=0; j<3; j++)
        {
            temp[j] = sampledU[i][j]; 
        }
        sampledValues[i] = temp;
    }

    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void Foam::SampledVelocityField::sample
(
    Foam::scalarListListList & sampledValues,
    const Foam::labelListList & indexList
) const
{
    Info<< "Sampling velocity for patch " << patch().name() << nl;
    
    const volVectorField & UField = mesh().lookupObject<volVectorField>("U");
    const vectorField & Uwall = UField.boundaryField()[patch().index()];
    

    forAll(indexList, i)
    {
        sampledValues[i] = scalarListList(indexList[i].size());
        forAll(indexList[i], j)
        {
            scalarList temp(3, 0.0);

            for (int k=0; k<3; k++)
            {
                temp[k] = 
                    (UField[indexList[i][j]] - Uwall[i])[k]; 
            }
            sampledValues[i][j] = temp;
        }
    }

    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void Foam::SampledVelocityField::registerFields
(
    const labelList &  indexList
) const
{
    // Initialize to 0
    scalarListList sampledU(patch().size());
    forAll(sampledU, i)
    {
        sampledU[i] = scalarList(3, 0.0);
    }
        
    if (mesh().foundObject<volVectorField>("U"))
    {
        const volVectorField & U = mesh().lookupObject<volVectorField>("U");

        forAll(sampledU, i)
        {
            forAll(sampledU[i], j)
            {
                sampledU[i][j] = U[indexList[i]][j];
            }
        }

        Helpers::projectOnPatch(patch().nf(), sampledU);
    }

    if (!db().foundObject<scalarListIOList>("U"))
    {
        mesh().time().store
        (        
            new IOList<scalarList>
            (
                IOobject
                (
                    "U",
                    mesh().time().timeName(),
                    db(), 
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                sampledU
            )
        );
    }
}


void Foam::SampledVelocityField::registerFields
(
    const labelListList & indexList
) const
{
    scalarListListList sampledU(patch().size());
    forAll(sampledU, i)
    {
        sampledU[i] = scalarListList(indexList[i].size());

        forAll(sampledU[i], j)
        {
            sampledU[i][j] = scalarList(3, 0.0);
        }
    }

    if (mesh().foundObject<volVectorField>("U"))
    {
        const auto & U = mesh().lookupObject<volVectorField>("U");
        forAll(sampledU, i)
        {
            forAll(sampledU[i], j)
            {
                forAll(sampledU[i][j], k)
                {
                    sampledU[i][j][k] = U[indexList[i][j]][k];
                }
            }
        }
        Helpers::projectOnPatch(patch().nf(), sampledU);
    }


    if (!db().foundObject<scalarListListIOList>("U"))
    {
        mesh().time().store
        (        
            new scalarListListIOList
            (
                IOobject
                (
                    "U",
                    mesh().time().timeName(),
                    db(), 
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                sampledU
            )
        );
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
