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

#include "SampledPGradField.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "List.H"
#include "helpers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SampledPGradField, 0);
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SampledPGradField::sample
(
    Foam::scalarListList & sampledValues,
    const Foam::labelList & indexList,
    const Foam::scalarField & h
) const
{
    if (debug)
    {
        Info<< "Sampling pressure gradient for patch " << patch_.name() << nl;
    }

    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField & faceNormals = tfaceNormals();
    
    const volVectorField & pGradField =
        mesh().lookupObject<volVectorField>("pGrad");
    autoPtr<interpolation<vector> > interpolator;
    interpolator.operator=
    (
        interpolation<vector>::New(interpolationType(), pGradField)
    );

    vectorField sampledPGrad(indexList.size());
    
    for (int i=0; i<indexList.size(); i++)
    {
        point p = faceCentres[i] - h[i]*faceNormals[i];
        const vector interp = interpolator->interpolate(p, indexList[i]);
        sampledPGrad[i] = interp;
        scalarList temp(3, 0.0);
        
        for (int j=0; j<3; j++)
        {
            temp[j] = sampledPGrad[i][j]; 
        }
        sampledValues[i] = temp;
    }
    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void
Foam::SampledPGradField::sample
(
    Foam::scalarListListList & sampledValues,
    const Foam::labelListList & indexList
) const
{
    Info<< "Sampling pressure gradient for patch " << patch().name() << nl;
    
    const volVectorField & pGradField =
        mesh().lookupObject<volVectorField>("pGrad");
    
    forAll(indexList, i)
    {
        sampledValues[i] = scalarListList(indexList[i].size());
        forAll(indexList[i], j)
        {
            sampledValues[i][j] = scalarList(3);

            forAll(sampledValues[i][j], k)
            {
                sampledValues[i][j][k] = pGradField[indexList[i][j]][k]; 
            }
        }
    }

    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void Foam::SampledPGradField::registerFields
(
    const labelList &  indexList

) const
{

    // Init sampled p grad to (0 0 0)
    scalarListList sampledPGrad(patch().size(), scalarList(3, 0.0));

    // Copy values from the pGrad field
    if (mesh().foundObject<volVectorField>("pGrad"))
    {
        forAll(sampledPGrad, i)
        {
            const volVectorField & pGrad =
                mesh().lookupObject<volVectorField>("pGrad");
            forAll(sampledPGrad[i], j)
            {
                sampledPGrad[i][j] = pGrad[indexList[i]][j];
            }
        }

        Helpers::projectOnPatch(patch().nf(), sampledPGrad);
    }

    
    if (!db().foundObject<scalarListIOList>("pGrad"))
    {
        mesh().thisDb().store
        (          
            new IOList<scalarList>
            (
                IOobject
                (
                    "pGrad",
                    mesh().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                sampledPGrad
             )
        );
    }
    
}


void Foam::SampledPGradField::registerFields
(
    const labelListList & indexList
) const
{
    scalarListListList sampledPGrad(patch().size());

    forAll(sampledPGrad, i)
    {
        sampledPGrad[i] = scalarListList(indexList[i].size());

        forAll(sampledPGrad[i], j)
        {
            sampledPGrad[i][j] = scalarList(3, 0.0);
        }
    }

    if (mesh().foundObject<volVectorField>("pGrad"))
    {
        const auto & pGrad = mesh().lookupObject<volVectorField>("pGrad");

        forAll(sampledPGrad, i)
        {
            forAll(sampledPGrad[i], j)
            {
                forAll(sampledPGrad[i][j], k)
                {
                    sampledPGrad[i][j][k] = pGrad[indexList[i][j]][k];
                }
            }
        }
        Helpers::projectOnPatch(patch().nf(), sampledPGrad);
    }
    
    if (!db().foundObject<scalarListIOList>("pGrad"))
    {
        mesh().thisDb().store
        (          
            new scalarListListIOList
            (
                IOobject
                (
                    "pGrad",
                    mesh().time().timeName(),
                    db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                sampledPGrad
             )
        );
    }
    
}

void Foam::SampledPGradField::recompute() const
{
    volVectorField & pGrad = const_cast<volVectorField &>
    (
        mesh().lookupObject<volVectorField>("pGrad")
    );
    const volScalarField & p = mesh().lookupObject<volScalarField>("p");
    
    pGrad = fvc::grad(p);
  
}


void Foam::SampledPGradField::createField
(
) const
{
    word hName;

    // Grab h for the current patch
    if (mesh_.foundObject<volScalarField>("hSampler"))
    {
        hName = "hSampler";
    }
    else
    {
        hName = "h";
    }

    const volScalarField & h = mesh_.lookupObject<volScalarField> (hName);
    
    if (!mesh().foundObject<volVectorField>("pGrad"))
    {
        mesh().thisDb().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "pGrad",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "pGrad",
                    dimLength/sqr(dimTime),
                    pTraits<vector>::zero
                ),
                h.boundaryField().types()
            )
        );
    }

    const IOdictionary fvSchemes = 
        mesh().thisDb().lookupObject<IOdictionary>("fvSchemes");

    // If we have p and scheme, compute the actual grad(p)
    // The check for the scheme is due to fvSchemes not read by decomposePar
    if ((mesh().foundObject<volScalarField>("p")) && 
        (fvSchemes.isDict("gradSchemes")))
    {
        recompute();
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
