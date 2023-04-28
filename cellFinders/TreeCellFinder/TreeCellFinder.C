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

    if (treePtr->nodes().empty() && (maxH != 0) && (patch().size()))
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
        // grab the wall-adjacent cell
        bool inside = 
#ifdef FOAM_VOLUMETYPE_NOT_CAPITAL
            boundaryTreePtr->getVolumeType(point) == volumeType::inside;
#else
            boundaryTreePtr->getVolumeType(point) == volumeType::INSIDE;
#endif

        if (h[i] < 0)
        {
            Warning
                << "TreeCellFinder: " << h[i]
                << " is negative and thus not a valid distance. "
                << "Will fall back to wall-adjacent cell for face "
                <<  i << " on patch " << patch().name() << nl;
                indexList[i] = faceCells[i];
        }
        else if ( (h[i] == 0) || (treePtr->nodes().empty()) )
        {
            indexList[i] = faceCells[i];
        }
        else if (!inside)
        {
            Warning
                << "TreeCellFinder: the point " << h[i]
                << " away from the wall is outside the domain. "
                << "Will fall back to wall-adjacent cell for face "
                <<  i << " on patch " << patch().name() << nl;
                indexList[i] = faceCells[i];
        }    
        else
        {

            pih = treePtr->findNearest(point, treePtr->bb().mag());
            indexList[i] = searchCellLabels[pih.index()];
        }
    }
    if (debug)
    {
        Info << "TreeCellFinder: Done" << nl;
    }
}

void Foam::TreeCellFinder::findCellIndices
(
    labelListList & indexList,
    const scalarField & h,
    const bool excludeWallAdjacent
) const
{
    scalar maxH = max(h);

    if (debug)
    {
        Info<< "TreeCellFinder: Constructing mesh bounding box" << nl;
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
    const volVectorField & C = mesh_.C();
    
    // Grab the global indices of adjacent cells 
    const UList<label> & faceCells = patch().faceCells();

    // Point 2h away from the face
    vector p;
    pointIndexHit pih;

    if (debug)
    {
        Info << "TreeCellFinder: Starting search for sampling cells" << nl;
    }

    forAll(faceCentres, i)
    {
        // Grab the point 2h away along the face normal
        // This just sets the direction of the ray cast
        p = faceCentres[i] - 2*faceNormals[i]*h[i];

        vector tolVector = (p - faceCentres[i])*1e-6;
        point startP = faceCentres[i] + tolVector;
        point endP = p + tolVector;

        // If h is zero, or the point is outside the domain,
        // set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly

        if (h[i] < 0)
        {
            Warning
                << "TreeCellFinder: " << h[i]
                << " is negative and thus not a valid distance. "
                << "Will fall back to wall-adjacent cell for face "
                <<  i << " on patch " << patch().name() << nl;
            indexList[i] = labelList(1, faceCells[i]);
        }
        else if ((h[i] == 0) || (treePtr->nodes().empty()))
        {
            indexList[i] = labelList(1, faceCells[i]);
        }
        else
        {
            // Allocate list for the sampling cell indices for face i
            indexList[i] = labelList(50);

            // Amount of cells sampled from for this face
            label n = 0;

            while (true)
            {
                if (debug > 1)
                {
                    Info << "Searching faces along line from: " << startP << " to " << endP << nl;
                }
                pih = treePtr->findLine(startP, endP);

                if (pih.hit())
                {
                    const point hitP = pih.hitPoint();
                    const label cellI =
                        treePtr->findNearest(hitP - tolVector, treePtr->bb().mag()).index();

                    indexList[i][n] = searchCellLabels[cellI];

                    if (debug > 1)
                    {
                        Info<< "Hit face: " << hitP << nl;
                        Info<< "CC: " << C[searchCellLabels[cellI]] << nl;
                    }
                   
                    n++;

                    // If we are now inside the last cell
                    if (mag(hitP - faceCentres[i]) >= h[i])
                    {
                        if (debug > 1)
                        {
                            Info<< "This is the last cell, n = " << n << nl;
                        }
                        indexList[i].setSize(n);
                        break;
                    }

                    startP = hitP + tolVector;
                }
                else
                {
                    if (n == 0)
                    {
                        // Not a single face intersected
                        // revert to wall-adjacent cell
                        Info << "No faces were intersected, reverting to "
                             << "wall-adjacent cell" << nl;
                        indexList[i].setSize(1, faceCells[i]);
                        break;
                    
                    }
                    else
                    {
                        // We went outside the domain
                        indexList[i].setSize(n);
                        if (debug > 1)
                        {
                            Info<< "This is the last cell, n = " << n << nl;
                        }
                        break;
                    }
                }
           }
        }

        // Remove wall adjacent cell from list if applicable
        if ((indexList[i].size() != 1) && (excludeWallAdjacent))
        {
            for(int j=0; j<indexList[i].size()-1; j++)
            {
                indexList[i][j] = indexList[i][j+1];
            }
            indexList[i].setSize(indexList[i].size()-1);
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

Foam::tmp<Foam::volScalarField> Foam::TreeCellFinder::distanceField() const
{
    

    word hName;

    // Grab h for the current patch
    if (mesh_.foundObject<volScalarField>("hSampler"))
    {
        hName = "hSampler";
    }
    else
    {
        hName = "h";
    }
    const volScalarField & h = mesh_.lookupObject<volScalarField>(hName);


    if (debug)
    {
        Info<< "CellFinder: Creating dist field" << nl;
    }

    tmp<volScalarField> dist
    ( 
        new volScalarField
        (
            IOobject
            (
                "dist",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dist", dimLength, 0),
            h.boundaryField().types()
        )
    );

    bool precomputedDist = 
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        mag(max(dist().primitiveField())) > VSMALL;
    #else
        mag(max(dist().internalField())) > VSMALL;
    #endif
    
    if (debug)
    {
        Info<<"CellFinder: using precumputed distance field" << nl;
    }

    if (!precomputedDist)
    {
        labelHashSet patchIDs(1);
        patchIDs.insert(patch().index());

        dictionary methodDict = dictionary();
        methodDict.lookupOrAddDefault(word("method"), word("meshWave"));

        if (debug)
        {
            Info<< "Initializing patchDistanceMethod" << nl;
        }

        autoPtr<patchDistMethod> pdm
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
            Info<< "CellFinder: Computing dist field" << nl;
        }
#ifdef FOAM_NEW_TMP_RULES
        pdm->correct(dist.ref());
#else
        pdm->correct(dist());
#endif
    }

    return dist;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
