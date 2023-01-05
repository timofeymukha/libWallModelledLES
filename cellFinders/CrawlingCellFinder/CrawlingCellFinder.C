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

#include "CrawlingCellFinder.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "codeRules.H"
#include <math.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CrawlingCellFinder, 0);
    addToRunTimeSelectionTable(CellFinder, CrawlingCellFinder, Patch);
}
#endif


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CrawlingCellFinder::CrawlingCellFinder
(
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "CrawlingCellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::CrawlingCellFinder::CrawlingCellFinder
(
    const word & cellFinderName,
    const fvPatch & p
)
:
   CellFinder(p)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::CrawlingCellFinder::findCellIndices
(
    labelList & indexList,
    const scalarField & h,
    const bool hIsIndex
) const
{
    const labelList & owner = mesh_.faceOwner();
    const labelList & neighbour = mesh_.faceNeighbour();
    const vectorField & faceCentres = mesh().Cf().primitiveField();
    const List<cell> & cells = mesh().cells();

    const vectorField & patchFaceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    const volVectorField & C = mesh_.C();
    const UList<label> & faceCells = patch().faceCells();
    const List<face> & faces = mesh().faces();
    const label patchStart = patch().start();
    const polyBoundaryMesh & boundaryMesh = mesh().boundaryMesh();

    if (hIsIndex)
    {
        Info<< "CrawlingCellFinder: Treating h as indices!" << nl;
    }
    else
    {
        Info<< "CrawlingCellFinder: Treating h as distances!" << nl;
    }

    forAll(patch(), patchFaceI)
    {

        label nLayers = 0;

        if (hIsIndex)
        {
            nLayers = std::round(h[patchFaceI]);
            if (nLayers < 1)
            {
                Warning
                    << "CrawlingCellFinder: Could not round "
                    << h[patchFaceI] << " to a valid cell layer index. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI] = faceCells[patchFaceI];
                continue;
            }
            else if (nLayers == 1)
            {
                indexList[patchFaceI] = faceCells[patchFaceI];
                continue;
            }
        }
        else
        {
            if (h[patchFaceI] < 0)
            {
                Warning
                    << "CrawlingCellFinder: " << h[patchFaceI]
                    << " is negative and thus not a valid distance. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI] = faceCells[patchFaceI];
                continue;
            }
            else if (h[patchFaceI] == 0)
            {
                indexList[patchFaceI] = faceCells[patchFaceI];
                continue;
            }

            nLayers = 100; // arbitrary big number
        }


        label startCellIndex = faceCells[patchFaceI];
        cell startCell = cells[startCellIndex];
        label startFaceLabel = patchStart + patchFaceI;
        label nextCellIndex = -1;

        for (label layer=0; layer < nLayers-1; layer++)
        {
            label opposingFace = 
                startCell.opposingFaceLabel(startFaceLabel, faces);

            if (opposingFace == -1)
            {
                Warning
                    << "CrawlingCellFinder: Could not find opposing face"
                    << "for cell " << startCellIndex << " with cell center "
                    << C[startCellIndex] << " and face " << startFaceLabel << nl
                    << "Will stop crawling and use the last valid cell " 
                    << "corresponding to index " << layer + 1
                    << nl;
                
                indexList[patchFaceI] = startCellIndex;
                break;
            }
            else if (opposingFace > mesh().nInternalFaces())
            {
                label opposingPatchInd = boundaryMesh.whichPatch(opposingFace);
                const fvPatch & opposingPatch =
                    mesh().boundary()[opposingPatchInd];
                word opposingPatchName = opposingPatch.name();

                // distance that will be used due to abortion in crawling
                scalar distance =
                    mag(C[startCellIndex] - patchFaceCentres[patchFaceI]);

                // If distance-based we might actually get the right h
                // in the last cell, so need to check for that 
                // before issuing the warning
                if (hIsIndex || mag(distance - h[patchFaceI])/distance > SMALL)
                {
                    Warning
                        << "CrawlingCellFinder: The opposing face for cell "
                        << startCellIndex << " with cell center "
                        << C[startCellIndex] << " and face " << startFaceLabel
                        << " belongs to patch " << opposingPatchName << nl 
                        << "Will stop crawling and use the last valid cell " 
                        << "corresponding to index " << layer + 1
                        << " and cell centre at distance " << distance
                        << ". Requested index/distance is " << h[patchFaceI]
                        << " which may or may not correspond to this cell." 
                        << nl;
                }

                break;
            }

            scalar distance =
                mag(faceCentres[opposingFace] - patchFaceCentres[patchFaceI]);

            // if the opposing face is above h, stop and grab the cell
            if ((!hIsIndex) && (distance > h[patchFaceI]))
            {
                break; 
            }

            if (owner[opposingFace] == startCellIndex)
            {
                nextCellIndex = neighbour[opposingFace];
            }
            else
            {
                nextCellIndex = owner[opposingFace];
            }

            startFaceLabel = opposingFace;
            startCellIndex = nextCellIndex;
            startCell = cells[startCellIndex];
        }

        indexList[patchFaceI] = startCellIndex;
    }
}


void Foam::CrawlingCellFinder::findCellIndices
(
    labelListList & indexList,
    const scalarField & h,
    const bool hIsIndex,
    const bool excludeWallAdjacent
) const
{
    const labelList & owner = mesh_.faceOwner();
    const labelList & neighbour = mesh_.faceNeighbour();
    const vectorField & faceCentres = mesh().Cf().primitiveField();
    const List<cell> & cells = mesh().cells();

    const vectorField & patchFaceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    const volVectorField & C = mesh_.C();
    const UList<label> & faceCells = patch().faceCells();
    const List<face> & faces = mesh().faces();
    const label patchStart = patch().start();
    const polyBoundaryMesh & boundaryMesh = mesh().boundaryMesh();

    if (hIsIndex)
    {
        Info<< "CrawlingCellFinder: Treating h as indices!" << nl;
    }
    else
    {
        Info<< "CrawlingCellFinder: Treating h as distances!" << nl;
    }

    forAll(patch(), patchFaceI)
    {

        label nLayers = 0;

        if (hIsIndex)
        {
            nLayers = std::round(h[patchFaceI]);
            if (nLayers < 1)
            {
                Warning
                    << "CrawlingCellFinder: Could not round "
                    << h[patchFaceI] << " to a valid cell layer index. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
            else if (nLayers == 1)
            {
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
        }
        else
        {
            if (h[patchFaceI] < 0)
            {
                Warning
                    << "CrawlingCellFinder: " << h[patchFaceI]
                    << " is negative and thus not a valid distance. "
                    << "Will fall back to wall-adjacent cell for face "
                    <<  patchFaceI << " on patch " << patch().name() << nl;
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }
            else if (h[patchFaceI] == 0)
            {
                indexList[patchFaceI].setSize(1, faceCells[patchFaceI]);
                continue;
            }

            nLayers = 100; // arbitrary big number
        }
        indexList[patchFaceI].setSize(nLayers, -1);

        label startCellIndex = faceCells[patchFaceI];
        cell startCell = cells[startCellIndex];
        label startFaceLabel = patchStart + patchFaceI;
        label nextCellIndex = -1;

        label layerCounter = 0;

        for (label layer=0; layer < nLayers; layer++)
        {
            layerCounter++;

            // Grab the next cell
            indexList[patchFaceI][layer] = startCellIndex;

            label opposingFace = 
                startCell.opposingFaceLabel(startFaceLabel, faces);

            if (opposingFace == -1)
            {
                Warning
                    << "CrawlingCellFinder: Could not find opposing face"
                    << "for cell " << startCellIndex << " with cell center "
                    << C[startCellIndex] << " and face " << startFaceLabel << nl
                    << "Will stop crawling and use the last valid cell " 
                    << "corresponding to index " << layer + 1
                    << nl;
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }
            // if we hit a patch
            else if (opposingFace > mesh().nInternalFaces())
            {
                label opposingPatchInd = boundaryMesh.whichPatch(opposingFace);
                const fvPatch & opposingPatch =
                    mesh().boundary()[opposingPatchInd];
                word opposingPatchName = opposingPatch.name();

                // Recompute as the actual distance that will be used
                // due to abortion in crawling
                scalar distance =
                    mag(C[startCellIndex] - patchFaceCentres[patchFaceI]);

                Warning
                    << "CrawlingCellFinder: The opposing face for cell "
                    << startCellIndex << " with cell center "
                    << C[startCellIndex] << " and face " << startFaceLabel
                    << " belongs to patch " << opposingPatchName << nl 
                    << "Will stop crawling and use the last valid cell " 
                    << "corresponding to index " << layer + 1
                    << " and distance " << distance
                    << nl;
                    
                // Need this when hIsIndex
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }

            scalar distance =
                mag(faceCentres[opposingFace] - patchFaceCentres[patchFaceI]);

            // if the opposing face is above h, stop 
            if ((!hIsIndex) && (distance > h[patchFaceI]))
            {
                indexList[patchFaceI].setSize(layerCounter);
                break;
            }

            // find index of the next cell
            if (owner[opposingFace] == startCellIndex)
            {
                nextCellIndex = neighbour[opposingFace];
            }
            else
            {
                nextCellIndex = owner[opposingFace];
            }

            startFaceLabel = opposingFace;
            startCellIndex = nextCellIndex;
            startCell = cells[startCellIndex];
        }

        if ((indexList[patchFaceI].size() != 1) && excludeWallAdjacent)
        {
            for(int j=0; j<indexList[patchFaceI].size()-1; j++)
            {
                indexList[patchFaceI][j] = indexList[patchFaceI][j+1];
            }
            indexList[patchFaceI].setSize(indexList[patchFaceI].size()-1);

        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
