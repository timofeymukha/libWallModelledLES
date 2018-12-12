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
#include "treeDataFace.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "codeRules.H"
#include "patchDistMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Sampler, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Sampler::createFields()
{
      
    if (!mesh_.thisDb().found("h"))
    {
        if (debug)
        {
            Info<< "Sampler: Creating h field" << nl;
        }

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

    // Try to read field with distance from cell to patches
    
    if (debug)
    {
        Info<< "Sampler: Creating dist field" << nl;
    }

    if (!mesh_.thisDb().found("dist"))
    {
        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "dist",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("dist", dimLength,0),
                h.boundaryField().types()
            )
        );
    }
    volScalarField & dist = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField>("dist"));

    if (mag(max(dist.internalField()).value()) < VSMALL)
    {
        if (debug)
        {
            Info<< "Sampler: Computing dist field" << nl;
        }

        labelHashSet patchIDs(1);
        patchIDs.insert(patch().index());
        //Info<<patchIDs << nl;

        dictionary methodDict = dictionary();
        methodDict.lookupOrAddDefault(word("method"), word("meshWave"));

        if (debug)
        {
            Info<< "    Initializing patchDistanceMethod" << nl;
        }

        autoPtr<patchDistMethod>  pdm_
        (
            patchDistMethod::New
            (
                methodDict,
                mesh_,
                patchIDs
            )
        );

        if (debug)
        {
            Info<< "    Computing" << nl;
        }
        pdm_->correct(dist);
    }
    
   // Field that marks cells that are used for sampling
    if (!mesh_.thisDb().found("samplingCells"))
    {
        if (debug)
        {
            Info<< "Sampler: Creating samplingCells field" << nl;
        }

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
    if (debug)
    {
        Info << "Sampler: The maximum h value is " << maxH << nl;
    }

    const volVectorField & C = mesh_.C();

    if (debug)
    {
        Info<< "Sampler: Constructing mesh bounding box" << nl;
    }

    treeBoundBox boundBox(mesh_.bounds());

    labelList searchCellLabels(C.size());
    label nSearchCells = 0;

    const volScalarField & distanceField = 
        mesh_.lookupObject<volScalarField> ("dist");

    if (debug)
    {
        Info<< "Sampler: Searching for cells closer to 2maxH to the wall"
            << nl;
    }

    forAll(searchCellLabels, i)
    {
        if (distanceField[i] < 2*maxH)
        {
            searchCellLabels[nSearchCells] = i;
            nSearchCells++;
        }
    }

    searchCellLabels.resize(nSearchCells);

    if (debug)
    {
        Info<< "Sampler: Found " << searchCellLabels.size() << " cells" << nl;
        Info<< "Sampler: Constructing cell octree" << nl;
    }
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

    if (treePtr->nodes().empty() && (maxH != 0))
    {
        Warning
            << "Sampler: max(h) is " << maxH << " but no cell centres within "
            << "distance 2*max(h) were found. "
            << "Will sample from wall-adjacent cells." << nl; 
    }


    if (debug)
    {
        Info<< "Sampler: Constructing face octree" << nl;
    }

    List<label> bndFaces(mesh_.nFaces() - mesh_.nInternalFaces());
    forAll(bndFaces, i)
    {
        bndFaces[i] = mesh_.nInternalFaces() + i;
    }

    indexedOctree<treeDataFace> * boundaryTreePtr =
    new indexedOctree<treeDataFace>
    (
        treeDataFace
        (
            false,
            mesh_,
            bndFaces
        ),
        boundBox,
        8,
        10,
        3.0
    );


    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the global indices of adjacent cells 
    const UList<label> & faceCells = patch().faceCells();

    vector point;
    //label pih;
    pointIndexHit pih;

    if (debug)
    {
        Info << "Sampler: Starting search for sampling cells" << nl;
    }
    
    forAll(faceCentres, i)
    {
        // Grab the point h away along the face normal
        point = faceCentres[i] - faceNormals[i]*h_[i];

        // If h is zero, or the point is outside the domain,
        // set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly
        bool inside = boundaryTreePtr->getVolumeType(point) == volumeType::INSIDE;
        if ((h_[i] == 0) || (!inside) || (treePtr->nodes().empty()))
        {
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList_[i] = faceCells[i];          
        }
        else
        {

            pih = treePtr->findNearest(point, treePtr->bb().mag());
            indexList_[i] = searchCellLabels[pih.index()];
            h_[i] = mag(C[indexList_[i]] - faceCentres[i]);
        }
    }
    if (debug)
    {
        Info << "Sampler: Done" << nl;
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
