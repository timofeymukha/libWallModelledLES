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

    scalarField hPatch = h.boundaryField()[patchIndex];

    scalar maxH = max(hPatch);

    if (debug)
    {
        Info<< "MultiCellSampler: Constructing mesh bounding box" << nl;
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
    const vectorField & faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField & cellCentres = tcellCentres();
    const volVectorField & C = mesh_.C();
    
    // Grab the global indices of adjacent cells 
    const UList<label> & faceCells = patch().faceCells();

    vector p;
    pointIndexHit pih;

    if (debug)
    {
        Info << "MultiCellSampler: Starting search for sampling cells" << nl;
    }
    
    forAll(faceCentres, i)
    {
        // Grab the point h away along the face normal
        p = faceCentres[i] - faceNormals[i]*hPatch[i];

        vector tolVector = (p - faceCentres[i])*1e-6;
        point startP = faceCentres[i] + tolVector;
        point endP = p + tolVector;

        // If h is zero, or the point is outside the domain,
        // set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly
        bool inside = boundaryTreePtr->getVolumeType(p) == volumeType::INSIDE;
        if ((hPatch[i] == 0) || (!inside) || (treePtr->nodes().empty()))
        {
            indexList_[i] = labelList(1);
            indexList_[i][0] = faceCells[i];

            h_[i] = scalarList(1);
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
        }
        else
        {
            // Allocate list for the sampling cell indices for face i, and associated heights
            indexList_[i] = labelList(50);
            h_[i] = scalarList(50);

            // Amount of cells sampled from for this face
            label n = 0;

            while (true)
            {
                Info << "Start: " << startP << " End: " << endP << nl;
                pih = treePtr->findLine(startP, endP);

                if (pih.hit())
                {
                    const point hitP = pih.hitPoint();
                    const label cellI =
                        treePtr->findNearest(hitP - tolVector, treePtr->bb().mag()).index();

                    indexList_[i][n] = searchCellLabels[cellI];
                    h_[i][n] = mag(C[searchCellLabels[cellI]] - faceCentres[i]);

                    Info<< "Hit face: " << hitP << nl;
                    Info<< "CC: " << C[searchCellLabels[cellI]] << nl;
                    

                    n++;

                    if (mag(hitP - faceCentres[i]) >= hPatch[i])
                    {
                        indexList_[i].setSize(n);
                        h_[i].setSize(n);
                        break;
                    }

                    startP = hitP + tolVector;
                }
                else
                {
                    // Not a single face intersected, revert to wall-adjacent cell
                    if (n == 0)
                    {
                        Info << "No faces were intersected, reverting to wall-adjacent cell" << nl;
                        indexList_[i].setSize(1);
                        indexList_[i][0] = faceCells[i];

                        h_[i].setSize(1);
                        h_[i] = mag(cellCentres[i] - faceCentres[i]);
                        break;
                    
                    }
                    indexList_[i].setSize(n);
                    h_[i].setSize(n);
                    break;
                }
            }

        }
        hPatch[i] = h_[i][h_[i].size() -1];
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
        hPatch;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            mesh_.lookupObject<volScalarField> ("samplingCells")
        );
    

    label totalSize = 0;
    forAll(indexList_, i)
    {
        totalSize += indexList_[i].size();

        forAll(indexList_[i], j)
        {
            samplingCells[indexList_[i][j]] = patchIndex;
        }
    }
    
    //TODO parallel
    label totalPatchSize =  patch().size();
    reduce(totalPatchSize, sumOp<label>());
    reduce(totalSize, sumOp<scalar>());
    if (totalPatchSize > 0)
    {
        Info<< "Average number of sampling cells per face is " <<
                totalSize/totalPatchSize << nl;
    }
}

void Foam::MultiCellSampler::createLengthList()
{
    // Cell volumes
    const scalarField & V = mesh_.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = scalarList(indexList_[i].size());
        
        forAll(lengthList_[i], j)
        {
            lengthList_[i][j] = pow(V[indexList_[i][j]], 1.0/3.0);
        }
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
    h_(p.size()),
    lengthList_(p.size())
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
        sampledFields_[fieldI]->sample(sampledList, indexList());
        
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
    field->registerFields(indexList());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
