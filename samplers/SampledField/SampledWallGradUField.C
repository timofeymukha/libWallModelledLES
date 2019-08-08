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

#include "SampledWallGradUField.H"
#include "volFields.H"
#include "codeRules.H"
#include "helpers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SampledWallGradUField, 0);
    addToRunTimeSelectionTable(SampledField, SampledWallGradUField, FvPatch);
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::SampledWallGradUField::sample
(
    Foam::scalarListList & sampledValues,
    const Foam::labelList & indexList
) const
{
    Info<< "Sampling wall-normal velocity gradient for patch "
        << patch_.name() << nl;
    
    label pI = patch().index();
   
    const volVectorField & wallGradU =
        mesh().lookupObject<volVectorField>("wallGradU");
    
    const vectorField & boundaryValues = wallGradU.boundaryField()[pI];
    
    for (int i=0; i<sampledValues.size(); i++)
    {
        scalarList temp(3, 0.0);
        
        for (int j=0; j<3; j++)
        {
            temp[j] = boundaryValues[i][j]; 
        }
        sampledValues[i] = temp;
    }

    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void
Foam::SampledWallGradUField::sample
(
    Foam::scalarListListList & sampledValues,
    const Foam::labelListList & indexListList
) const
{
    Info<< "Sampling wall-normal velocity gradient for patch "
        << patch_.name() << nl;
    
    label pI = patch().index();
   
    const volVectorField & wallGradU =
        mesh().lookupObject<volVectorField>("wallGradU");
    
    const vectorField & boundaryValues = wallGradU.boundaryField()[pI];
    
    for (int i=0; i<indexListList.size(); i++)
    {
        sampledValues[i] = scalarListList(1);
        sampledValues[i][0] = scalarList(3);
        
        for (int j=0; j<3; j++)
        {
            sampledValues[i][0][j] = boundaryValues[i][j]; 
        }
    }

    Helpers::projectOnPatch(patch().nf(), sampledValues);
}


void Foam::SampledWallGradUField::registerFields
(
    const labelList &  indexList 
) const
{
    // Grab h to copy bcs from it.
    const volScalarField & h = mesh().lookupObject<volScalarField>("h");
    
    if (!mesh().foundObject<volVectorField>("wallGradU"))
    {
        mesh().time().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "wallGradU",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "wallGradU",
                    dimVelocity/dimLength,
                    pTraits<vector>::zero
                ),
                h.boundaryField().types()
            )
        );  
    }
    
    scalarListList sampledWallGradU(patch().size());

    forAll(sampledWallGradU, i)
    {
        sampledWallGradU[i] = scalarList(3, 0.0);
    }

    if (mesh().foundObject<volVectorField>("U"))
    {
        recompute();
        
        const volVectorField & wallGradU = 
            mesh().lookupObject<volVectorField>("wallGradU");
        
        label pI = patch().index();
        const vectorField & boundaryValues = wallGradU.boundaryField()[pI];

        forAll(sampledWallGradU, i)
        {
            forAll(sampledWallGradU[i], j)
            {
                sampledWallGradU[i][j] = boundaryValues[i][j];
            }
        }

        Helpers::projectOnPatch(patch().nf(), sampledWallGradU);
    }

    
    mesh().time().store
    (        
        new IOList<scalarList>
        (
            IOobject
            (
                "wallGradU",
                mesh().time().timeName(),
                db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sampledWallGradU
        )
    );

}

void Foam::SampledWallGradUField::registerFields
(
    const labelListList & indexListList
) const
{
    // Grab h to copy bcs from it.
    const volScalarField & h = mesh().lookupObject<volScalarField>("h");
    
    if (!mesh().foundObject<volVectorField>("wallGradU"))
    {
        mesh().time().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "wallGradU",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "wallGradU",
                    dimVelocity/dimLength,
                    pTraits<vector>::zero
                ),
                h.boundaryField().types()
            )
        );  
    }

    const objectRegistry & registry =
        mesh().subRegistry("wallModelSampling").subRegistry(patch().name());

    scalarListListList sampledWallGradU(patch().size());

    if (mesh().foundObject<volVectorField>("U"))
    {
        
        forAll(sampledWallGradU, i)
        {
            sampledWallGradU[i] = scalarListList(1);
            sampledWallGradU[i][0] = scalarList(3, 0.0);
        }
    }

    
    mesh().time().store
    (        
        new scalarListListIOList
        (
            IOobject
            (
                "wallGradU",
                mesh().time().timeName(),
                registry,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sampledWallGradU
        )
    );

}


void Foam::SampledWallGradUField::recompute() const
{
    label pI = patch().index(); 
    
    volVectorField & wallGradU = const_cast<volVectorField &>
    (
            mesh().lookupObject<volVectorField>("wallGradU")
    );
    
    const volVectorField & U = mesh().lookupObject<volVectorField>("U");
    const fvPatchVectorField & Uwall = U.boundaryField()[pI];
      
    vectorField Udiff(Uwall.patchInternalField() - Uwall);
#ifdef FOAM_NEW_GEOMFIELD_RULES
    wallGradU.boundaryFieldRef()[pI]
#else        
    wallGradU.boundaryField()[pI]
#endif
    ==
        patch().deltaCoeffs()*Udiff;  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
