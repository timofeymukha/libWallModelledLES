#include "fvCFD.H"
#include "CellFinder.H"
#include "CrawlingCellFinder.H"
#include "NormalRayCellFinder.H"
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

    const wordList names({"CrawlingCellFinder", "NormalRayCellFinder"});

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

    // The removed Tree sampler points to the NormalRay sampler
    ASSERT_DEATH
    (
        {
            CellFinder::New("TreeCellFinder", patch);
        },
        "NormalRay"
    );
}


TEST_F(CellFinderTest, FindCellsSerial)
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
    NormalRayCellFinder normalRayFinder(patch);

    auto check = [&](const CellFinder & finder)
    {
        List<CellFinder::samplingLocation> locations;
        mapDistributeBase map;
        finder.findCells(locations, map, h, false);

        // In serial, all locations are on this processor
        ASSERT_EQ(map.constructSize(), patch.size());
        ASSERT_EQ(locations.size(), patch.size());

        labelList cells(locations.size());

        forAll(locations, i)
        {
            const label faceI = locations[i].originFace;

            ASSERT_EQ(locations[i].originProc, Pstream::myProcNo());
            ASSERT_EQ(locations[i].wallCentre, patch.Cf()[faceI]);
            ASSERT_FLOAT_EQ(locations[i].h, 0.45);
            ASSERT_FLOAT_EQ(mesh.C()[locations[i].cell][1], 0.5);
            cells[i] = locations[i].cell;
        }

        // The map sends the values to the corresponding faces
        map.distribute(cells);

        forAll(cells, faceI)
        {
            ASSERT_FLOAT_EQ(mesh.C()[cells[faceI]][1], 0.5);
            ASSERT_FLOAT_EQ
            (
                mesh.C()[cells[faceI]][0],
                patch.Cf()[faceI][0]
            );
        }
    };

    check(crawlingFinder);
    check(normalRayFinder);
}


TEST_F(CellFinderTest, NormalRayHIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    autoPtr<CellFinder> finder(CellFinder::New("NormalRayCellFinder", patch));

    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;

    ASSERT_DEATH
    (
        {
            finder->findCells
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
