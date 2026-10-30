#include "fvCFD.H"
#include "CrawlingParallelCellFinder.H"
#include "CrawlingCellFinder.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class CrawlingParallelCellFinderTest : public ChannelFlow
{};


// Collect the found cells in the order of the patch faces
static labelList crawlingCellsByFace
(
    const List<CellFinder::samplingLocation> & locations,
    const mapDistributeBase & map
)
{
    labelList cells(locations.size());

    forAll(locations, i)
    {
        cells[i] = locations[i].cell;
    }

    map.distribute(cells);
    return cells;
}


TEST_F(CrawlingParallelCellFinderTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingParallelCellFinder finder(patch);

    ASSERT_EQ(&finder.patch(), &patch);
    ASSERT_EQ(&finder.mesh(), &mesh);
}


TEST_F(CrawlingParallelCellFinderTest, Locations)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    const vectorField & faceCentres = patch.Cf();
    const tmp<vectorField> tfaceNormals = patch.nf();

    CrawlingParallelCellFinder finder(patch);

    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;

    finder.findCells
    (
        locations,
        map,
        scalarField(patch.size(), 3),
        true
    );

    ASSERT_EQ(locations.size(), patch.size());
    ASSERT_EQ(map.constructSize(), patch.size());

    boolList foundFace(patch.size(), false);

    forAll(locations, i)
    {
        const CellFinder::samplingLocation & location = locations[i];
        const label faceI = location.originFace;

        ASSERT_EQ(location.originProc, Pstream::myProcNo());
        ASSERT_EQ(location.layer, 3);
        ASSERT_EQ(location.wallCentre, faceCentres[faceI]);
        ASSERT_EQ(location.wallNormal, tfaceNormals()[faceI]);
        ASSERT_FLOAT_EQ(location.h, 3);
        ASSERT_FALSE(foundFace[faceI]);
        foundFace[faceI] = true;
    }
}


TEST_F(CrawlingParallelCellFinderTest, IndexMatchesCrawlingCellFinder)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const volVectorField & C = mesh.C();

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch & patch = mesh.boundary()[patchI];

        if (patch.name() != "bottomWall" && patch.name() != "topWall")
        {
            continue;
        }

        CrawlingParallelCellFinder finder(patch);
        CrawlingCellFinder crawlingFinder(patch);

        // Includes invalid and too large indices
        const scalarList indices({-1, 0, 1, 2, 5, 6, 10, 15});

        forAll(indices, indexI)
        {
            const scalarField h(patch.size(), indices[indexI]);

            List<CellFinder::samplingLocation> locations;
            mapDistributeBase map;
            finder.findCells(locations, map, h, true);

            labelList expected(patch.size(), -1);
            crawlingFinder.findCellIndices(expected, h, true);

            const labelList cells(crawlingCellsByFace(locations, map));

            forAll(cells, i)
            {
                ASSERT_EQ(cells[i], expected[i]);
            }
        }
    }

    // Sanity check on the expected cells
    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingParallelCellFinder finder(patch);
    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;
    finder.findCells(locations, map, scalarField(patch.size(), 2), true);

    const labelList cells(crawlingCellsByFace(locations, map));

    forAll(cells, i)
    {
        ASSERT_FLOAT_EQ(C[cells[i]][1], 0.3);
    }
}


TEST_F(CrawlingParallelCellFinderTest, DistanceMatchesCrawlingCellFinder)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const volVectorField & C = mesh.C();

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch & patch = mesh.boundary()[patchI];

        if (patch.name() != "bottomWall" && patch.name() != "topWall")
        {
            continue;
        }

        CrawlingParallelCellFinder finder(patch);
        CrawlingCellFinder crawlingFinder(patch);

        // Includes invalid and too large distances
        const scalarList distances({-0.1, 0, 0.1, 0.25, 0.5, 1.05, 1.9, 5});

        forAll(distances, distanceI)
        {
            const scalarField h(patch.size(), distances[distanceI]);

            List<CellFinder::samplingLocation> locations;
            mapDistributeBase map;
            finder.findCells(locations, map, h, false);

            labelList expected(patch.size(), -1);
            crawlingFinder.findCellIndices(expected, h, false);

            const labelList cells(crawlingCellsByFace(locations, map));

            forAll(cells, i)
            {
                ASSERT_EQ(cells[i], expected[i]);
            }
        }
    }

    // Sanity check on the expected cells
    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingParallelCellFinder finder(patch);
    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;
    finder.findCells
    (
        locations,
        map,
        scalarField(patch.size(), 1.05),
        false
    );

    const labelList cells(crawlingCellsByFace(locations, map));

    forAll(cells, i)
    {
        ASSERT_FLOAT_EQ(C[cells[i]][1], 1.1);
    }
}
