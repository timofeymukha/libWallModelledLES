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

#include "Sampler.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "IOField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::Sampler::createFields() const
{
   
    //const objectRegistry & db = patch_.boundaryMesh().mesh();
    
    if (!mesh_.found("h"))
    {
        mesh_.store
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch_.boundaryMesh().mesh()
            )
        );
    }
    
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));
    
   // Field that marks cells that are used for sampling
    if (!mesh_.found("samplingCells"))
    {
        mesh_.store
        (
            new volScalarField
            (
                IOobject
                (
                    "samplingCells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("samplingCells", dimless,0),
                h.boundaryField().types()
            )
        );
    }
    
    // Field for sampled velocity
    new IOField<vector>
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            db_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField(patch_.size(), pTraits<vector>::zero)
    );

}


void Foam::Sampler::createIndexList()
{
   
    const label patchIndex = patch().index();
    
    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh.lookupObject<volScalarField> ("h"));
    
    h_ = h.boundaryField()[patchIndex];



    // Create a searcher for the mesh
    meshSearch ms(mesh);
    
    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the global indices of adjacent cells 
    const labelUList & faceCells = patch().faceCells();

    vector point;
    forAll(faceCentres, i)
    {
        // If h is zero, set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly.
        if (h_[i] == 0)
        {          
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList_[i] = faceCells[i];          
        }
        else
        {
            // Grab the point h away along the face normal
            point = faceCentres[i] - faceNormals[i]*h_[i];

            // Check that point is inside the (processor) domain
            // Otherwise fall back to adjacent cell's centre.
            if (!ms.isInside(point))
            {
                point = cellCentres[i];
            }

            // Find the cell where the point is located
            indexList_[i] = ms.findNearestCell(point, -1, true);
                        
            // Set h to the distance between face centre and located cell's
            // center
            h_[i] = mag(mesh.C()[indexList_[i]] - faceCentres[i]);
        }
    }
    
    // Assign computed h_ to the global h field
    h.boundaryField()[patch().index()] == h_;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            mesh.lookupObject<volScalarField> ("samplingCells")
        );
    
    forAll(indexList_, i)
    {
        samplingCells[indexList_[i]] = patchIndex; 
    }
}


void Foam::Sampler::createLengthList()
{
    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh(); 
    
    // Cell volumes
    const scalarField & V = mesh.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = pow(V[indexList_[i]], 1.0/3.0);
    }
    
}

void Foam::Sampler::project(vectorField & field) const
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Sampler::Sampler
(
    const fvPatch & p,
    scalar averagingTime
)
:
    patch_(p),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0),
    averagingTime_(averagingTime),
    db_(patch_.boundaryMesh().mesh().subRegistry("wallModelSampling", 1)),
    mesh_(patch_.boundaryMesh().mesh())
{
    createFields();
    createIndexList();
    createLengthList();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Sampler::sample() const
{

    // Weight for time-averaging, default to 1 i.e no averaging.
    scalar eps = 1;
    if (averagingTime_ > mesh_.time().deltaTValue())
    {
        eps = mesh_.time().deltaTValue()/averagingTime_;
    }
    
    
    forAll(db_.names(), fieldNameI)
    {
        word name = db_.names()[fieldNameI];
        
        // Sample velocity
        const volVectorField & field = mesh_.lookupObject<volVectorField>(name);
        const vectorField & internal = field.internalField();
        //const fvPatchVectorField & wall = Ufield.boundaryField()[patch_.index()];

        vectorField & sampled = 
                const_cast<vectorField &>(db_.lookupObject<vectorField>(name));
        
        vectorField Up(patch().size()); 

        // NO SUPPORT FOR MOVING WALLS FOR NOW
        forAll(Up, i)
        {   
        //    Up[i] = Uinternal[indexList_[i]] - Uwall[i];
              Up[i] = internal[indexList_[i]];
        }

        project(Up);

        forAll(sampled, i)  
        {    
            sampled[i] = eps*Up[i] + (1 - eps)*sampled[i];
        }
    }
    
    
    
    
    
    
    
    /*
     * 
     * 
     * 
     * 
            // Sample velocity
        const volVectorField & Ufield = mesh_.lookupObject<volVectorField>("U");
        const vectorField & Uinternal = Ufield.internalField();
        const fvPatchVectorField & Uwall = Ufield.boundaryField()[patch_.index()];

        vectorField & Usampled = 
                const_cast<vectorField &>(db_.lookupObject<vectorField>("U"));
        vectorField Up(patch().size()); 

        forAll(Up, i)
        {   
            Up[i] = Uinternal[indexList_[i]] - Uwall[i];
        }

        project(Up);

        forAll(Usampled, i)  
        {    
            Usampled[i] = eps*Up[i] + (1 - eps)*Usampled[i];
        }
    const volScalarField & p = db().lookupObject<volScalarField>("p");

    //calculate pressure Gradient
    vectorField gradPField = fvc::grad(p);
    vectorField gradP(patch().size());
    
    forAll (gradP, i)
    {
        gradP[i] = gradPField[cellIndexList_[i]];
    }
    
    project(gradP);
    
    
    
    forAll (pressureGrad_, i)
    {
        pressureGrad_[i] = eps*gradP[i] + (1 - eps)*pressureGrad_[i];
    }
      
    // Wall-normal velocity gradient
    vectorField Udiff = Uwall.patchInternalField() - Uwall;
    project(Udiff);
    const vectorField wallGradU(patch().deltaCoeffs()*Udiff);
    
    volVectorField & wallGradUField = 
        const_cast<volVectorField &>
        (
            db().lookupObject<volVectorField>("wallGradU")
        );
    wallGradUField.boundaryField()[patch().index()] == wallGradU;

    forAll(U_, i)  
    {    
        wallGradU_[i] = eps*wallGradU[i] + (1 - eps)*wallGradU_[i];
    }
    */

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //