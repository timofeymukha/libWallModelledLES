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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarListList Foam::SampledWallGradUField::sample() const
{
    
    Info << "Sampling wall-normal velocity gradient" << nl;
    
    label pI = patch().index();
    scalarListList sampledValues(cellIndexList_.size());
    
    const volVectorField & wallGradU =
        db().lookupObject<volVectorField>("wallGradU");
    
    const vectorField & boundaryValues = wallGradU.boundaryField()[pI];
    
    for(int i=0; i<cellIndexList_.size(); i++)
    {
        sampledValues[i] = List<scalar>(3);
        
        for(int j=0; j<3; j++)
        {
            sampledValues[i][j] = boundaryValues[i][j]; 
        }
    }
    projectVectors(sampledValues);
    return sampledValues;
}

void Foam::SampledWallGradUField::registerFields() const
{
    // Grab h to copy bcs from it.
    const volScalarField & h = db().lookupObject<volScalarField>("h");
    
    db().thisDb().store
    (     
        new volVectorField
        (
            IOobject
            (
                "wallGradU",
                db().time().timeName(),
                db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            db(),
            dimensionedVector
            (
                "wallGradU",
                dimVelocity/dimLength,
                vector(0, 0, 0)
            ),
            h.boundaryField().types()
        )
    );  
    
    db().thisDb().store
    (        
        new IOField<vector>
        (
            IOobject
            (
                "wallGradU",
                db().time().timeName(),
                db().subRegistry("wallModelSampling", 0),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            vectorField(cellIndexList_.size(), pTraits<vector>::zero)
        )
    );
}

void Foam::SampledWallGradUField::recompute() const
{
    label pI = patch().index(); 
    
    volVectorField & wallGradU = const_cast<volVectorField &>
    (
            db().lookupObject<volVectorField>("wallGradU")
    );
    
    const volVectorField & U = db().lookupObject<volVectorField>("U");
    const fvPatchVectorField & Uwall = U.boundaryField()[pI];
      
    vectorField Udiff = Uwall.patchInternalField() - Uwall;
    //project(Udiff);

    wallGradU.boundaryField()[pI] == patch().deltaCoeffs()*Udiff;  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //