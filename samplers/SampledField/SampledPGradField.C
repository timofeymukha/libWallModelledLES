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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SampledPGradField::sample
(
    Foam::scalarListList & sampledValues,
    const Foam::labelList & indexList
) const
{
    Info<< "Sampling pressure gradient for patch " << patch_.name() << nl;
    
    const volVectorField & pGradField =
        db().lookupObject<volVectorField>("pGrad");
    vectorField sampledPGrad(indexList.size());
    
    for (int i=0; i<indexList.size(); i++)
    {
        sampledPGrad[i] = pGradField[indexList[i]];
        sampledValues[i] = List<scalar>(3);
        
        for (int j=0; j<3; j++)
        {
            sampledValues[i][j] = sampledPGrad[i][j]; 
        }
    }
    projectVectors(sampledValues);
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
        db().lookupObject<volVectorField>("pGrad");
    
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
    projectVectors(sampledValues);
}


void Foam::SampledPGradField::registerFields() const
{
    // Grab h to copy bcs from it.
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    if (!db().foundObject<volVectorField>("pGrad"))
    {
        db().thisDb().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "pGrad",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                db(),
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
    
    const objectRegistry & registry =
        db().subRegistry("wallModelSampling").subRegistry(patch_.name());

    scalarListList sampledPGrad(patch().size());

    if (db().thisDb().foundObject<volScalarField>("p"))
    {
        //recompute();
        
        forAll(sampledPGrad, i)
        {
            sampledPGrad[i] = scalarList(3, 0.0);
            //forAll(sampledPGrad[i], j)
            //{
                //sampledPGrad[i][j] = 0; //pGrad[cellIndexList_[i]][j];
            //}
        }
        
        //projectVectors(sampledPGrad);
    }
    
    db().thisDb().store
    (          
        new IOList<scalarList>
        (
            IOobject
            (
                "pGrad",
                db().time().timeName(),
                registry,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sampledPGrad
         )
    );
    
}


void Foam::SampledPGradField::registerFields
(
    const labelListList & indexListList
) const
{
    // Grab h to copy bcs from it.
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    if (!db().foundObject<volVectorField>("pGrad"))
    {
        db().thisDb().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "pGrad",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                db(),
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

    const objectRegistry & registry =
        db().subRegistry("wallModelSampling").subRegistry(patch_.name());

    scalarListListList sampledPGrad(patch().size());

    if (db().thisDb().foundObject<volScalarField>("p"))
    {
        forAll(sampledPGrad, i)
        {
            sampledPGrad[i] = scalarListList(indexListList[i].size());

            forAll(sampledPGrad[i], j)
            {
                sampledPGrad[i][j] = scalarList(3, 0.0);
            }
        }
    }
    
    db().thisDb().store
    (          
        new scalarListListIOList
        (
            IOobject
            (
                "pGrad",
                db().time().timeName(),
                registry,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sampledPGrad
         )
    );
    
}

void Foam::SampledPGradField::recompute() const
{
    volVectorField & pGrad = const_cast<volVectorField &>
    (
            db().lookupObject<volVectorField>("pGrad")
    );
    const volScalarField & p = db().lookupObject<volScalarField>("p");
    
    pGrad = fvc::grad(p);
  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
