#include "fvCFD.H"
#include "CellFinder.H"
#include "CrawlingCellFinder.H"
#include "TreeCellFinder.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class CellFinderTest : public ChannelFlow
{};


TEST_F(CellFinderTest, RunTimeSelection)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    const wordList names
    ({
        "CrawlingCellFinder",
        "TreeCellFinder",
        "CrawlingParallelCellFinder",
        "TreeParallelCellFinder"
    });

    forAll(names, i)
    {
        autoPtr<CellFinder> finder(CellFinder::New(names[i], patch));
        ASSERT_EQ(finder->type(), names[i]);
        ASSERT_EQ(&finder->patch(), &patch);
    }

    ASSERT_DEATH
    (
        {
            CellFinder::New("RandomNameCellFinder", patch);
        },
        "FATAL ERROR"
    );
}


TEST_F(CellFinderTest, LocalFindCells)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    const scalarField h(patch.size(), 0.45);

    CrawlingCellFinder crawlingFinder(patch);
    TreeCellFinder treeFinder(patch);

    labelList crawlingCells(patch.size(), -1);
    crawlingFinder.findCellIndices(crawlingCells, h, false);

    labelList treeCells(patch.size(), -1);
    treeFinder.findCellIndices(treeCells, h);

    auto check = [&](const CellFinder & finder, const labelList & expected)
    {
        List<CellFinder::samplingLocation> locations;
        mapDistributeBase map;
        finder.findCells(locations, map, h, false);

        // The map does not involve other processors
        ASSERT_EQ(map.comm(), UPstream::commSelf());
        ASSERT_EQ(map.constructSize(), patch.size());
        ASSERT_EQ(locations.size(), patch.size());

        forAll(locations, i)
        {
            ASSERT_EQ(locations[i].originProc, Pstream::myProcNo());
            ASSERT_EQ(locations[i].originFace, i);
            ASSERT_EQ(locations[i].cell, expected[i]);
            ASSERT_EQ(locations[i].wallCentre, patch.Cf()[i]);
            ASSERT_FLOAT_EQ(locations[i].h, 0.45);
        }
    };

    check(crawlingFinder, crawlingCells);
    check(treeFinder, treeCells);
}


TEST_F(CellFinderTest, TreeHIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    TreeCellFinder finder(patch);

    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;

    ASSERT_DEATH
    (
        {
            finder.findCells
            (
                locations,
                map,
                scalarField(patch.size(), 2),
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(CellFinderTest, SamplingLocationIO)
{
    CellFinder::samplingLocation location;
    location.originProc = 1;
    location.originFace = 2;
    location.cell = 3;
    location.face = 4;
    location.layer = 5;
    location.wallCentre = point(1, 2, 3);
    location.wallNormal = vector(0, -1, 0);
    location.h = 0.5;

    List<CellFinder::samplingLocation> locations(2, location);
    locations[1].cell = 7;

    OStringStream os;
    os << locations;

    IStringStream is(os.str());
    List<CellFinder::samplingLocation> readLocations(is);

    ASSERT_EQ(readLocations.size(), 2);
    ASSERT_TRUE(readLocations[0] == location);
    ASSERT_TRUE(readLocations[1] == locations[1]);
    ASSERT_TRUE(readLocations[0] != readLocations[1]);
}
