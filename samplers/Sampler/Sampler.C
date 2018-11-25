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
#include "objectRegistry.H"
#include "IOField.H"
#include "SampledField.H"
#include "SampledVelocityField.H"
#include "SampledPGradField.H"
#include "SampledWallGradUField.H"
#include "treeDataCell.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "codeRules.H"
#include "wallDist.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Sampler::createFields()
{
      
    if (!mesh_.thisDb().found("h"))
    {
        mesh_.thisDb().store
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
                mesh_
            )
        );
    }
    
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));
    
   // Field that marks cells that are used for sampling
    if (!mesh_.thisDb().found("samplingCells"))
    {
        mesh_.thisDb().store
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
      
}


void Foam::Sampler::createIndexList()
{
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));
    
    h_ = h.boundaryField()[patchIndex];

    scalar maxH = max(h_);

    const volVectorField & C = mesh_.C();

    treeBoundBox boundBox(mesh_.bounds());


    wallDist distance(mesh_);
    const volScalarField & distanceField = distance.y();


    // Indices of cells lying closer than 2h to the patch
    labelList searchCellLabels(C.size());
    label nSearchCells = 0;
    forAll(searchCellLabels, i)
    {
        if (distanceField[i] < 2*maxH)
        {
            searchCellLabels[nSearchCells] = i;
            nSearchCells++;
        }
    }
    searchCellLabels.resize(nSearchCells);

    Info << "Planting the tree" << nl;
    indexedOctree<treeDataCell> * treePtr = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,
            mesh_,
            searchCellLabels,
            polyMesh::CELL_TETS
        ),
        boundBox,
        8,
        10,
        3.0
    );
    Info << "Grown" << nl;

    // Create a searcher for the mesh
    Info << "Creating meshsearch" << nl;
    meshSearch ms(mesh_);
    Info << "Done" << nl;
    
    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the global indices of adjacent cells 
    const UList<label> & faceCells = patch().faceCells();

    vector point;
    pointIndexHit pih;

    forAll(faceCentres, i)
    {
        // If h is zero, set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly
        // Note that if maxH is 0, we will never use our empty tree
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
//            Info << "Check if inside" << nl;
//            if (!ms.isInside(point))
//            {
//                point = cellCentres[i];
//            }
//            Info << "Done" << nl;

            // Find the cell where the point is located
//            indexList_[i] = ms.findNearestCell(point, -1, true);
//            pih = treePtr->findNearest(point, treePtr->bb().mag());
            pih = treePtr->findNearest(point, treePtr->bb().mag());
            indexList_[i] = searchCellLabels[pih.index()];
//            Info << indexList_[i] << " " << pih << " " << searchCellLabels[pih.index()] << nl;
            Info << treePtr->findInside(point) << nl;
                        
            // Set h to the distance between face centre and located cell's
            // center
            h_[i] = mag(C[indexList_[i]] - faceCentres[i]);
        }
    }
    
    // Assign computed h_ to the global h field
#ifdef FOAM_NEW_GEOMFIELD_RULES
    h.boundaryFieldRef()[patch().index()]
#else        
    h.boundaryField()[patch().index()]
#endif
    ==
        h_;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            mesh_.lookupObject<volScalarField> ("samplingCells")
        );
    
    forAll(indexList_, i)
    {
        samplingCells[indexList_[i]] = patchIndex; 
    }
}


void Foam::Sampler::createLengthList()
{
    // Cell volumes
    const scalarField & V = mesh_.V();
    
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
    mesh_(patch_.boundaryMesh().mesh()),
    sampledFields_(0)
{
    if (!mesh_.foundObject<objectRegistry>("wallModelSampling"))
    {
        objectRegistry * subObr = new objectRegistry
        (
            IOobject
            (
                "wallModelSampling",
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        subObr->store();
    }

    objectRegistry * subObr = new objectRegistry
    (
        IOobject
        (
            patch_.name(),
            mesh_.time().constant(),
            mesh_.subRegistry("wallModelSampling"),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    subObr->store();

    createFields();
    createIndexList();
    createLengthList();
    
    addField
    (
            new SampledVelocityField(patch_, indexList_)     
    );
    
    addField
    (
            new SampledWallGradUField(patch_, indexList_)     
    );
}

Foam::Sampler::Sampler(const Sampler & copy)
:
    patch_(copy.patch_),
    indexList_(copy.indexList_),
    lengthList_(copy.lengthList_),
    h_(copy.h_),
    averagingTime_(copy.averagingTime_),
    mesh_(copy.mesh_),
    sampledFields_(copy.sampledFields_.size())
{
    forAll(copy.sampledFields_, i)
    {
        sampledFields_[i] = copy.sampledFields_[i]->clone();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Sampler::~Sampler()
{
    forAll(sampledFields_, i)
    {
        delete sampledFields_[i];
    }
}


void Foam::Sampler::sample() const
{
    // Ensure this processor has part of the patch
    if (!patch().size())
    {
        return;
    }    
    
    // Weight for time-averaging, default to 1 i.e no averaging.
    scalar eps = 1;
    if (averagingTime_ > mesh_.time().deltaTValue())
    {
        eps = mesh_.time().deltaTValue()/averagingTime_;
    }

    forAll(sampledFields_, fieldI)
    {

        scalarListList sampledList(patch().size());
        sampledFields_[fieldI]->sample(sampledList);
        
        if (sampledFields_[fieldI]->nDims() == 3)
        {
            vectorField sampledField(patch_.size());
            listListToField<vector>(sampledList, sampledField);
            
            vectorField & storedValues = const_cast<vectorField &>
            (
                db().lookupObject<vectorField>(sampledFields_[fieldI]->name())
            );
            storedValues = eps*sampledField + (1 - eps)*storedValues;
        }     

    }
}


template<class Type>
void Foam::Sampler::listListToField
(
    const scalarListList & list,
    Field<Type> & field    
) const
{
    forAll(list, i)
    {
        Type element;
        forAll(list[i], j)
        {
            element[j] = list[i][j];
        }
        field[i] = element;
    }
}

void Foam::Sampler::addField(SampledField * field)
{
//    sampledFields_.append(field);
    sampledFields_.setSize(sampledFields_.size() + 1, field);
}


void Foam::Sampler::recomputeFields() const
{  
    forAll(sampledFields_, i)
    {
        sampledFields_[i]->recompute();
    } 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
