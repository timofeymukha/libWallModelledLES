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

#include "MultiCellSampler.H"
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
#include "scalarListListIOList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MultiCellSampler, 0);
    addToRunTimeSelectionTable(Sampler, MultiCellSampler, PatchAndAveragingTime);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::MultiCellSampler::createIndexList()
{
    const label patchIndex = patch().index();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));

    
    h_ = h.boundaryField()[patchIndex];
    scalar maxH = max(h_);

    if (debug)
    {
        Info<< "MultiCellSampler: Constructing mesh bounding box" << nl;
    }

    treeBoundBox boundBox(mesh_.bounds());

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
            << "MultiCellSampler: max(h) is " << maxH << " but no cell centres within "
            << "distance 2*max(h) were found. "
            << "Will sample from wall-adjacent cells." << nl; 
    }


    if (debug)
    {
        Info<< "MultiCellSampler: Constructing face octree" << nl;
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

    vector p;
    //label pih;
    pointIndexHit pih;

    if (debug)
    {
        Info << "MultiCellSampler: Starting search for sampling cells" << nl;
    }
    
    forAll(faceCentres, i)
    {
        // Grab the point h away along the face normal
        p = faceCentres[i] - faceNormals[i]*h_[i];

        vector tolVector = (p - faceCentres[i])*1e-6;
        point startP = faceCentres[i] - tolVector;
        point endP = p + tolVector;

        // If h is zero, or the point is outside the domain,
        // set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly
        bool inside = boundaryTreePtr->getVolumeType(p) == volumeType::INSIDE;
        if ((h_[i] == 0) || (!inside) || (treePtr->nodes().empty()))
        {
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList_[i] = faceCells[i];          
            multiindexList_[i] = labelList(1);
            multiindexList_[i][0] = faceCells[i];
        }
        else
        {
            // Allocate list for the sampling cell indices for face i
            multiindexList_[i] = List<label>(50);

            // Amount of cells sampled from for this face
            label n = 0;

            while (true)
            {
                pih = treePtr->findLine(startP, endP);

                if (pih.hit())
                {
                    const point hitP = pih.hitPoint();
                    const label cellI =
                        treePtr->findNearest(hitP - tolVector, treePtr->bb().mag()).index();

                    multiindexList_[i][n] = searchCellLabels[cellI];
                    Info<< "Face: " << hitP << nl;
                    Info<< "CC: " << C[searchCellLabels[cellI]] << nl;
                    

                    n++;

                    if (mag(hitP - faceCentres[i]) >= h_[i])
                    {
                        Info<< mag(hitP - faceCentres[i]) << "  " << nl;
                        multiindexList_[i].setSize(n);
                        break;
                    }

                    startP = hitP + tolVector;
                    Info << startP << endP << nl;
                }
                else
                {
                    multiindexList_[i].setSize(n);
                    break;
                }
            }

            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList_[i] = faceCells[i];          
        }
    }
    if (debug)
    {
        Info << "MultiCellSampler: Done" << nl;
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
        //samplingCells[indexList_[i]] = patchIndex; 
    }

    label totalSize = 0;
    forAll(multiindexList_, i)
    {
        totalSize += multiindexList_[i].size();

        forAll(multiindexList_[i], j)
        {
            samplingCells[multiindexList_[i][j]] = patchIndex;
        }
    }
    
    //TODO parallel
    Info<< "Average number of sampling cells per face is " <<
        totalSize/patch().size() << nl;
}

void Foam::MultiCellSampler::createLengthList()
{
    // Cell volumes
    const scalarField & V = mesh_.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = pow(V[indexList_[i]], 1.0/3.0);
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MultiCellSampler::MultiCellSampler
(
    const fvPatch & p,
    scalar averagingTime
)
:
    Sampler(p, averagingTime),
    indexList_(p.size()),
    multiindexList_(p.size()),
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

    //List<vectorList> sampledU(10);
    //forAll (sampledU, i)
    //{
        //sampledU[i] = List<vector>(3);

        //forAll (sampledU[i], j)
        //{
            //sampledU[i][j] = pTraits<vector>::zero;
        //}
    //}
    
    //mesh().thisDb().store
    //(
        //new IOList<vectorList >
        //(
            //IOobject
            //(
                //"Utest",
                //db().time().timeName(),
                //mesh(), 
                //IOobject::READ_IF_PRESENT,
                //IOobject::AUTO_WRITE
            //),
            //sampledU
        //)
    //);

}


Foam::MultiCellSampler::MultiCellSampler
(
    const word & samplerName,
    const fvPatch & p,
    scalar averagingTime
)
:
    MultiCellSampler(p, averagingTime)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::MultiCellSampler::~MultiCellSampler()
{
}


void Foam::MultiCellSampler::sample() const
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

    Info << "Sampling" << nl;
    forAll(sampledFields_, fieldI)
    {

        scalarListListList sampledList(patch().size());
        sampledFields_[fieldI]->sample(sampledList, multiindexList());
        
        scalarListListIOList & storedValues = const_cast<scalarListListIOList & >
        (
            db().lookupObject<scalarListListIOList>(sampledFields_[fieldI]->name())
        );

        forAll(storedValues, i)
        {
            forAll(storedValues[i], j)
            {
                forAll(storedValues[i][j], k)
                storedValues[i][j][k] = eps*sampledList[i][j][k] +
                                     (1 - eps)*storedValues[i][j][k];
            }
        }

    }
}


void Foam::MultiCellSampler::addField(SampledField * field)
{
    Sampler::addField(field);
    field->registerFields(multiindexList());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
