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
#include "codeRules.H"
#include "patchDistMethod.H"
#include "scalarListIOList.H"
#include "scalarListListIOList.H"
#include "TreeCellFinder.H"
#include "CrawlingCellFinder.H"
#include "Sampler.H"
#include "surfaceMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(MultiCellSampler, 0);
    addToRunTimeSelectionTable(Sampler, MultiCellSampler, SamplerRTSTable);
}
#endif

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::MultiCellSampler::createIndexList()
{
    const label patchIndex = patch().index();

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
    
    volScalarField & hField = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> (hName));

    scalarField hPatch = hField.boundaryField()[patchIndex];

    if (cellFinderType() == "Crawling")
    {
        CrawlingCellFinder cellFinder(patch());
        cellFinder.findCellIndices(indexList_, hPatch, hIsIndex_, excludeWallAdjacent_);
    }
    else if (cellFinderType() == "Tree")
    {
        if (hIsIndex_)
        {
            FatalErrorInFunction
                << "MultiCellSampler: hIsIndex is not supported by the Tree "
                << "sampler. Please use the Crawling sampler or provide h as "
                << "distances."
                <<  abort(FatalError);
        }

        TreeCellFinder cellFinder(patch());
        cellFinder.findCellIndices(indexList_, hPatch, excludeWallAdjacent_);
    }
    else
    {
        FatalErrorInFunction
            << "MultiCellSampler: invalid sampler finder name, choose "
            << "Tree or Crawling. Current choice is " << cellFinderType()
            <<  abort(FatalError);
    }


    const vectorField & patchFaceCentres = patch().Cf();
    const volVectorField & C = mesh_.C();

    forAll(patch(), faceI)
    {
        h_[faceI].setSize(indexList_[faceI].size());
        forAll(indexList_[faceI], i)
        {
            h_[faceI][i] =
                mag(C[indexList_[faceI][i]] - patchFaceCentres[faceI]);
        }
    }

    scalarField hTop(indexList_.size());
    forAll(patch(), faceI)
    {
        const label n = h_[faceI].size() - 1;
        hTop[faceI] = h_[faceI][n];
    }


    // If the h field holds the distance, reassign the real distance used
    if (!hIsIndex())
    {
#ifdef FOAM_NEW_GEOMFIELD_RULES
        hField.boundaryFieldRef()[patch().index()]
#else
        hField.boundaryField()[patch().index()]
#endif
        ==
            hTop;
    }
    
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

void Foam::MultiCellSampler::createLengthList(const word lengthScaleType)
{
    if (lengthScaleType == "CubeRootVol")
    {
        createLengthListCubeRootVol();
    }
    else if (lengthScaleType == "WallNormalDistance")
    {
        createLengthListWallNormalDistance(); 
    }
    
}

void Foam::MultiCellSampler::createLengthListCubeRootVol()
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

void Foam::MultiCellSampler::createLengthListWallNormalDistance()
{
    const vectorField & faceCentres = mesh().Cf().primitiveField();
    const List<cell> & cells = mesh().cells();
    const vectorField & patchFaceCentres = patch().Cf();
    const List<face> & faces = mesh().faces();
    forAll(lengthList_, i)
    {
        lengthList_[i] = scalarList(indexList_[i].size());
        const vector patchFaceI = patchFaceCentres[i];

        forAll(lengthList_[i], j)
        {
            const label index = indexList_[i][j];
            const cell cellI = cells[index];

            scalar minDist = GREAT;
            label minDistFace = 0;

            for (int k=0; k<cellI.nFaces(); k++)
            {
                const vector faceK = faceCentres[cellI[k]];

                const scalar dist = mag(faceK - patchFaceI);
                if (dist < minDist)
                {
                    minDist = dist;
                    minDistFace = cellI[k];
                }
            }
            
            label opposingFace = 
                cellI.opposingFaceLabel(minDistFace, faces);
            vector opposingFaceCentre = faceCentres[opposingFace];
            vector minDistFaceCentre = faceCentres[minDistFace];

            lengthList_[i][j] = mag(opposingFaceCentre - minDistFaceCentre); 
        }

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MultiCellSampler::MultiCellSampler
(
    const fvPatch & p,
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
:
    Sampler(p, averagingTime, interpolationType, cellFinderType,
            lengthScaleType, hIsIndex, excludeWallAdjacent),
    indexList_(p.size()),
    h_(p.size()),
    lengthList_(p.size())
{

    if (interpolationType != "cell")
    {
        FatalErrorIn 
        (
            "MultiCellSampler::MultiCellSampler"
        )   << "MulticellSmapler: interpolation is not supported "
            << " for multicell sampling. Use 'cell'"
            << exit(FatalError); 
    }

    createIndexList();
    createLengthList(lengthScaleType);
    
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
    scalar averagingTime,
    const word interpolationType,
    const word cellFinderType,
    const word lengthScaleType,
    bool hIsIndex,
    bool excludeWallAdjacent
)
:
    MultiCellSampler
    (
        p,
        averagingTime,
        interpolationType,
        cellFinderType,
        lengthScaleType,
        hIsIndex,
        excludeWallAdjacent
    )
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    forAll(sampledFields_, fieldI)
    {

        scalarListListList sampledList(patch().size());
        sampledFields_[fieldI].sample(sampledList, indexList());
        
        scalarListListIOList & storedValues = const_cast<scalarListListIOList & >
        (
            db().lookupObject<scalarListListIOList>(sampledFields_[fieldI].name())
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
