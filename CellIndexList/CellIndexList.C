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

#include "CellIndexList.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::CellIndexList::createFields() const
{
   
    const objectRegistry & db = patch_.boundaryMesh().mesh();
    
    if (!db.found("h"))
    {
        db.store
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    db.time().timeName(),
                    db,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch_.boundaryMesh().mesh()
            )
        );
    }
    
    volScalarField & h = 
        const_cast<volScalarField &>(db.lookupObject<volScalarField> ("h"));
    
        // Field that marks cells that are used for sampling
    if (!db.found("samplingCells"))
    {
        db.store
        (
            new volScalarField
            (
                IOobject
                (
                    "samplingCells",
                    db.time().timeName(),
                    db,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("samplingCells", dimless,0),
                h.boundaryField().types()
            )
        );
    }
}


void Foam::CellIndexList::createIndexList()
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


void Foam::CellIndexList::createLengthList()
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CellIndexList::CellIndexList
(
    const fvPatch& p
)
:
    patch_(p),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0)
{
    createFields();
    createIndexList();
    createLengthList();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //