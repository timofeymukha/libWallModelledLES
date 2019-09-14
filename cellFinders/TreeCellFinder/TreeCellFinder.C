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

#include "TreeCellFinder.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "codeRules.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(TreeCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, TreeCellFinder, Patch);
}
#endif


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TreeCellFinder::TreeCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "TreeCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::TreeCellFinder::TreeCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::TreeCellFinder::findCellIndices
(
    labelList & indexList,
    const scalarField & h
) const
{
    const label patchIndex = patch().index();
    
    scalar maxH = max(h);

    if (debug)
    {
        Info<< "TreeCellSampler: Constructing mesh bounding box" << nl;
    }

    treeBoundBox boundBox(mesh().bounds());
    Random rndGen(261782);    
#ifdef FOAM_TREEBOUNDBOX_DOES_NOT_ACCEPT_RNG
    boundBox.extend(1e-4);
#else
    boundBox.extend(rndGen, 1e-4);
#endif
    boundBox.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    boundBox.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    tmp<Foam::volScalarField> distField(distanceField());

    const scalarField & dist =
#ifdef FOAM_NEW_GEOMFIELD_RULES
        distField().primitiveField();
#else
        distField().internalField();
#endif
    

    tmp<labelField> tSearchCellLabels = findCandidateCellLabels(dist, h);
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
            << "TreeCellFinder: max(h) is " << maxH << " but no cell centres within "
            << "distance 2*max(h) were found. "
            << "Will sample from wall-adjacent cells." << nl; 
    }


    if (debug)
    {
        Info<< "TreeCellFinder: Constructing face octree" << nl;
    }

    labelList bndFaces(mesh().nFaces() - mesh().nInternalFaces());
    forAll(bndFaces, i)
    {
        bndFaces[i] = mesh().nInternalFaces() + i;
    }

    autoPtr<indexedOctree<treeDataFace> > boundaryTreePtr
    (
        new indexedOctree<treeDataFace>
        (
            treeDataFace
            (
                false,
                mesh(),
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
    pointIndexHit pih;

    if (debug)
    {
        Info << "TreeCellFinder: Starting search for sampling cells" << nl;
    }
    
    forAll(faceCentres, i)
    {
        // Grab the point h away along the face normal
        point = faceCentres[i] - faceNormals[i]*h[i];

        // If h is zero, or the point is outside the domain,
        // set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly
        bool inside = 
#ifdef FOAM_VOLUMETYPE_NOT_CAPITAL
            boundaryTreePtr->getVolumeType(point) == volumeType::inside;
#else
            boundaryTreePtr->getVolumeType(point) == volumeType::INSIDE;
#endif
        if ((h[i] == 0) || (!inside) || (treePtr->nodes().empty()))
        {
            //h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList[i] = faceCells[i];
        }
        else
        {

            pih = treePtr->findNearest(point, treePtr->bb().mag());
            indexList[i] = searchCellLabels[pih.index()];
            //h_[i] = mag(C[indexList_[i]] - faceCentres[i]);
        }
    }
    if (debug)
    {
        Info << "TreeCellFinder: Done" << nl;
    }
}


Foam::tmp<Foam::labelField>
Foam::TreeCellFinder::findCandidateCellLabels
(
    const scalarField & dist,
    const scalarField & h
) const
{
    scalar maxH = max(h);
    if (debug)
    {
        Info << "TreeCellFinder: The maximum h value is " << maxH << nl;
    }

    if (debug)
    {
        Info<< "TreeCellFinderSampler: Constructing mesh bounding box" << nl;
    }

    tmp<labelField> tCandidates(new labelField(mesh().nCells()));
#ifdef FOAM_NEW_TMP_RULES
    labelField & candidates = tCandidates.ref();
#else
    labelField & candidates = tCandidates();
#endif

    label nCandidates = 0;

    if (debug)
    {
        Info<< "TreeCellFinder: Searching for cells closer to 2maxH to the wall"
            << nl;
    }

    forAll(candidates, i)
    {
        if (dist[i] < 2*maxH)
        {
            candidates[nCandidates] = i;
            nCandidates++;
        }
    }

    candidates.resize(nCandidates);

    if (debug)
    {
        Info<< "TreeCellFinder: Found " << candidates.size() << " candidate "
            << "cells" << nl;
    }

    return tCandidates;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
