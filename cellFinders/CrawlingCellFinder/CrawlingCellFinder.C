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
    const scalarField & h
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


    forAll(patch(), patchFaceI)
    {

        const label nLayers = std::round(h[patchFaceI]);

        if (nLayers <= 1)
        {
            indexList[patchFaceI] = faceCells[patchFaceI];
            break;
        }

        point patchFaceCenterI = patchFaceCentres[patchFaceI];
        Info << "patchFaceCenter " << patchFaceCenterI <<nl;

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
                    << C[startCellIndex] << "and face " << startFaceLabel << nl
                    << "Will stop crawling and use last valid cell for sampling"
                    << nl;
                
                indexList[patchFaceI] = startCellIndex;
                break;
            }
            else if (opposingFace > mesh().nInternalFaces())
            {
                Warning
                    << "CrawlingCellFinder: Opposing face for cell "
                    << startCellIndex << " with cell center "
                    << C[startCellIndex] << "and face " << startFaceLabel
                    << "belongs to a patch" << nl 
                    << "Will stop crawling and use last valid cell for sampling"
                    << nl;

                indexList[patchFaceI] = startCellIndex;
                break;
            }

            Info << "opposingFace " << opposingFace << faceCentres[opposingFace] << nl;

            if (owner[opposingFace] == startCellIndex)
            {
                nextCellIndex = neighbour[opposingFace];
            }
            else
            {
                nextCellIndex = owner[opposingFace];
            }

            Info << "nextCellIndex " << nextCellIndex << " " << C[nextCellIndex] << nl;
            startFaceLabel = opposingFace;
            startCellIndex = nextCellIndex;
            startCell = cells[startCellIndex];
        }

        indexList[patchFaceI] = startCellIndex;

        Info<< nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
