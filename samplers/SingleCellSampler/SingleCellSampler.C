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

#include "SingleCellSampler.H"
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
#include "scalarListIOList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SingleCellSampler, 0);
    addToRunTimeSelectionTable(Sampler, SingleCellSampler, PatchAndAveragingTime);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::SingleCellSampler::createIndexList()
{
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));

    
    h_ = h.boundaryField()[patchIndex];
    scalar maxH = max(h_);

    if (debug)
    {
        Info<< "SingleCellSampler: Constructing mesh bounding box" << nl;
    }

    treeBoundBox boundBox(mesh_.bounds());
    Random rndGen(261782);    
#ifdef FOAM_TREEBOUNDBOX_DOES_NOT_ACCEPT_RNG
    boundBox.extend(1e-4);
#else
    boundBox.extend(rndGen, 1e-4);
#endif
    boundBox.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    boundBox.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    tmp<labelField> tSearchCellLabels = findSearchCellLabels();
    const labelField & searchCellLabels = tSearchCellLabels();

    autoPtr<indexedOctree<treeDataCell> > treePtr
    (
        new indexedOctree<treeDataCell>
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
        )
    );

    if (treePtr->nodes().empty() && (maxH != 0))
    {
        Warning
            << "SingleCellSampler: max(h) is " << maxH << " but no cell centres within "
            << "distance 2*max(h) were found. "
            << "Will sample from wall-adjacent cells." << nl; 
    }


    if (debug)
    {
        Info<< "SingleCellSampler: Constructing face octree" << nl;
    }

    labelList bndFaces(mesh_.nFaces() - mesh_.nInternalFaces());
    forAll(bndFaces, i)
    {
        bndFaces[i] = mesh_.nInternalFaces() + i;
    }

    autoPtr<indexedOctree<treeDataFace> > boundaryTreePtr
    (
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
        )
    );


    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    const volVectorField & C = mesh_.C();
    
    // Grab the global indices of adjacent cells 
    const UList<label> & faceCells = patch().faceCells();

    vector point;
    //label pih;
    pointIndexHit pih;

    if (debug)
    {
        Info << "SingleCellSampler: Starting search for sampling cells" << nl;
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
        Info << "SingleCellSampler: Done" << nl;
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


void Foam::SingleCellSampler::createLengthList()
{
    // Cell volumes
    const scalarField & V = mesh_.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = pow(V[indexList_[i]], 1.0/3.0);
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SingleCellSampler::SingleCellSampler
(
    const fvPatch & p,
    scalar averagingTime
)
:
    Sampler(p, averagingTime),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0)
{
    createIndexList();
    createLengthList();
    
    addField
    (
            new SampledVelocityField(patch_)     
    );
    
    addField
    (
            new SampledWallGradUField(patch_)     
    );
}


Foam::SingleCellSampler::SingleCellSampler
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime
)
:
    SingleCellSampler(p, averagingTime)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::SingleCellSampler::~SingleCellSampler()
{
}


void Foam::SingleCellSampler::sample() const
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
        sampledFields_[fieldI]->sample(sampledList, indexList());
        
        scalarListIOList & storedValues = const_cast<scalarListIOList & >
        (
            db().lookupObject<scalarListIOList>(sampledFields_[fieldI]->name())
        );

        forAll(storedValues, i)
        {
            forAll(storedValues[i], j)
            {
                storedValues[i][j] = eps*sampledList[i][j] +
                                     (1 - eps)*storedValues[i][j];
            }
        }

    }
}


void Foam::SingleCellSampler::addField(SampledField * field)
{
    Sampler::addField(field);
    field->registerFields(indexList());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
