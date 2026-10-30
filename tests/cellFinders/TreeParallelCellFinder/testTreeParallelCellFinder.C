#include "fvCFD.H"
#include "TreeParallelCellFinder.H"
#include "TreeCellFinder.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class TreeParallelCellFinderTest : public ChannelFlow
{};


// Collect the found cells in the order of the patch faces
static labelList treeCellsByFace
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


TEST_F(TreeParallelCellFinderTest, MatchesTreeCellFinder)
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

        TreeParallelCellFinder finder(patch);
        TreeCellFinder treeFinder(patch);

        // Distances away from the cell faces, includes invalid and outside
        // the domain
        const scalarList distances({-0.1, 0, 0.1, 0.25, 0.5, 1.05, 1.9, 5});

        forAll(distances, distanceI)
        {
            const scalarField h(patch.size(), distances[distanceI]);

            List<CellFinder::samplingLocation> locations;
            mapDistributeBase map;
            finder.findCells(locations, map, h, false);

            labelList expected(patch.size(), -1);
            treeFinder.findCellIndices(expected, h);

            const labelList cells(treeCellsByFace(locations, map));

            forAll(cells, i)
            {
                ASSERT_EQ(cells[i], expected[i]);
            }
        }
    }

    // Sanity check on the expected cells
    const fvPatch & patch = mesh.boundary()["bottomWall"];
    TreeParallelCellFinder finder(patch);
    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;
    finder.findCells
    (
        locations,
        map,
        scalarField(patch.size(), 1.05),
        false
    );

    const labelList cells(treeCellsByFace(locations, map));

    forAll(cells, i)
    {
        ASSERT_FLOAT_EQ(C[cells[i]][1], 1.1);
    }
}


TEST_F(TreeParallelCellFinderTest, OutsideDomain)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    TreeParallelCellFinder finder(patch);

    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;
    finder.findCells(locations, map, scalarField(patch.size(), 5), false);

    const labelList cells(treeCellsByFace(locations, map));

    forAll(cells, i)
    {
        ASSERT_EQ(cells[i], patch.faceCells()[i]);
    }
}


TEST_F(TreeParallelCellFinderTest, HIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    TreeParallelCellFinder finder(patch);

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
