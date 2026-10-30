#include "fvCFD.H"
#include "NormalRayCellFinder.H"
#include "CrawlingCellFinder.H"
#include "pointIOField.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class NormalRayCellFinderTest : public ChannelFlow
{};


// Shear the mesh of the channel in the streamwise direction, x += a*y, which
// makes the cells skewed, while the walls remain flat
static void shearMesh(const Time & runTime, const scalar a)
{
    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(points, i)
    {
        points[i].x() += a*points[i].y();
    }

    points.write();
}


// Find the sampling locations in the order of the patch faces
static List<CellFinder::samplingLocation> locationsByFace
(
    const CellFinder & finder,
    const scalarField & h
)
{
    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;
    finder.findCells(locations, map, h, false);

    List<CellFinder::samplingLocation> byFace(locations.size());

    forAll(locations, i)
    {
        byFace[locations[i].originFace] = locations[i];
    }

    return byFace;
}


TEST_F(NormalRayCellFinderTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);

    ASSERT_EQ(&finder.patch(), &patch);
    ASSERT_EQ(&finder.mesh(), &mesh);
}


TEST_F(NormalRayCellFinderTest, FindDistanceBasedBottomWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    const scalarList hs({0.05, 0.1, 0.32, 1.9, 1.95});
    const scalarList expected({0.1, 0.1, 0.3, 1.9, 1.9});

    forAll(hs, hI)
    {
        findCellIndices
        (
            finder,
            indexList,
            scalarField(patch.size(), hs[hI]),
            false
        );

        forAll(indexList, i)
        {
            ASSERT_FLOAT_EQ(C[indexList[i]][1], expected[hI])
                << "h = " << hs[hI];
        }
    }
}


TEST_F(NormalRayCellFinderTest, FindDistanceBasedTopWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    NormalRayCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    const scalarList hs({0.05, 0.1, 0.32, 1.9, 1.95});
    const scalarList expected({1.9, 1.9, 1.7, 0.1, 0.1});

    forAll(hs, hI)
    {
        findCellIndices
        (
            finder,
            indexList,
            scalarField(patch.size(), hs[hI]),
            false
        );

        forAll(indexList, i)
        {
            ASSERT_FLOAT_EQ(C[indexList[i]][1], expected[hI])
                << "h = " << hs[hI];
        }
    }
}


// If the point is on the face between two cells, the cell further away from
// the wall is used, as by the CrawlingCellFinder
TEST_F(NormalRayCellFinderTest, FindDistanceOnFace)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    findCellIndices(finder, indexList, scalarField(patch.size(), 0.2), false);

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.3);
    }
}


TEST_F(NormalRayCellFinderTest, FindDistanceInvalidDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    findCellIndices(finder, indexList, scalarField(patch.size(), -0.1), false);

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(NormalRayCellFinderTest, FindDistanceZeroDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    findCellIndices(finder, indexList, scalarField(patch.size(), 0), false);

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


// The line leaves the domain through the opposite wall, the last cell is used
// and the sampling height is reduced to the middle of the line within it
TEST_F(NormalRayCellFinderTest, FindDistanceTooLargeDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);

    const List<CellFinder::samplingLocation> locations
    (
        locationsByFace(finder, scalarField(patch.size(), 10))
    );

    forAll(locations, i)
    {
        ASSERT_FLOAT_EQ(C[locations[i].cell][1], 1.9);
        ASSERT_FLOAT_EQ(locations[i].h, 1.9);
        ASSERT_EQ(locations[i].layer, 10);
    }
}


// The point is on the opposite wall, which is still within the domain
TEST_F(NormalRayCellFinderTest, FindDistanceOnOppositeWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);

    const List<CellFinder::samplingLocation> locations
    (
        locationsByFace(finder, scalarField(patch.size(), 2))
    );

    forAll(locations, i)
    {
        ASSERT_FLOAT_EQ(C[locations[i].cell][1], 1.9);
        ASSERT_FLOAT_EQ(locations[i].h, 2);
    }
}


TEST_F(NormalRayCellFinderTest, Locations)
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

    NormalRayCellFinder finder(patch);

    List<CellFinder::samplingLocation> locations;
    mapDistributeBase map;

    finder.findCells
    (
        locations,
        map,
        scalarField(patch.size(), 0.45),
        false
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
        ASSERT_FLOAT_EQ(mesh.C()[location.cell][1], 0.5);
        ASSERT_EQ(location.wallCentre, faceCentres[faceI]);
        ASSERT_EQ(location.wallNormal, tfaceNormals()[faceI]);
        ASSERT_FLOAT_EQ(location.h, 0.45);
        ASSERT_FALSE(foundFace[faceI]);
        foundFace[faceI] = true;
    }
}


TEST_F(NormalRayCellFinderTest, HIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);

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


// On a skewed mesh, the cell containing the point at distance h along the
// wall normal is found, whereas crawling follows the skewed columns of cells
TEST_F(NormalRayCellFinderTest, SkewedMeshCellContainsPoint)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    shearMesh(runTime, 2);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    label nFound = 0;
    label nOutside = 0;
    label nDifferentFromCrawling = 0;

    for (const word patchName : {"bottomWall", "topWall"})
    {
        const fvPatch & patch = mesh.boundary()[patchName];
        const tmp<vectorField> tfaceNormals = patch.nf();

        NormalRayCellFinder finder(patch);
        CrawlingCellFinder crawlingFinder(patch);

        for (const scalar h : {0.1, 0.35, 0.77, 1.3, 1.99})
        {
            const List<CellFinder::samplingLocation> locations
            (
                locationsByFace(finder, scalarField(patch.size(), h))
            );

            labelList crawlingCells;
            findCellIndices
            (
                crawlingFinder,
                crawlingCells,
                scalarField(patch.size(), h),
                false
            );

            forAll(locations, i)
            {
                const point p = patch.Cf()[i] - h*tfaceNormals()[i];
                const label cellI = mesh.findCell(p);

                if (cellI == -1)
                {
                    // The point is outside the domain, behind a cyclic patch
                    ASSERT_LT(locations[i].h, h);
                    nOutside++;
                }
                else
                {
                    ASSERT_EQ(locations[i].cell, cellI)
                        << patchName << " h = " << h << " face " << i;
                    ASSERT_FLOAT_EQ(locations[i].h, h);
                    nFound++;

                    if (crawlingCells[i] != cellI)
                    {
                        nDifferentFromCrawling++;
                    }
                }
            }
        }
    }

    ASSERT_GT(nFound, 0);
    ASSERT_GT(nOutside, 0);
    ASSERT_GT(nDifferentFromCrawling, 0);
}


TEST_F(NormalRayCellFinderTest, FindDistanceBasedBottomWallMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelListList indexList(patch.size());

    // Single, wall-adjacent cell
    finder.findCellIndices(indexList, scalarField(patch.size(), 0.1), false);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    // h is on the face between the first two cells, expecting a 2-cell list
    finder.findCellIndices(indexList, scalarField(patch.size(), 0.2), false);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 2);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
    }

    // h above 3rd cell center, expecting a 3-cell list
    finder.findCellIndices(indexList, scalarField(patch.size(), 0.55), false);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 3);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
        ASSERT_FLOAT_EQ(C[indexList[i][2]][1], 0.5);
    }
}


TEST_F(NormalRayCellFinderTest, FindDistanceBasedBottomWallWithExclusionMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelListList indexList(patch.size());

    // Single, wall-adjacent cell, which is kept
    finder.findCellIndices(indexList, scalarField(patch.size(), 0.1), true);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    // h above 3rd cell center, expecting a 2-cell list
    finder.findCellIndices(indexList, scalarField(patch.size(), 0.55), true);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 2);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.3);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.5);
    }
}


TEST_F(NormalRayCellFinderTest, ZeroAndInvalidDistanceMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelListList indexList(patch.size());

    for (const scalar h : {0.0, -10.0})
    {
        finder.findCellIndices(indexList, scalarField(patch.size(), h), true);

        forAll(indexList, i)
        {
            ASSERT_EQ(indexList[i].size(), 1);
            ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        }
    }
}


TEST_F(NormalRayCellFinderTest, TooLargeDistanceMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);
    labelListList indexList(patch.size());

    finder.findCellIndices(indexList, scalarField(patch.size(), 2.9), false);

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 10);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][9]][1], 1.9);
    }
}


// The multi-cell list ends with the single sampling cell, and consists of
// distinct cells, also on a skewed mesh
TEST_F(NormalRayCellFinderTest, SkewedMeshMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    shearMesh(runTime, 2);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    NormalRayCellFinder finder(patch);

    const scalarField h(patch.size(), 0.65);

    labelListList indexList(patch.size());
    finder.findCellIndices(indexList, h, false);

    labelList cells;
    findCellIndices(finder, cells, h, false);

    forAll(indexList, i)
    {
        ASSERT_GE(indexList[i].size(), 4);
        ASSERT_EQ(indexList[i].first(), patch.faceCells()[i]);
        ASSERT_EQ(indexList[i].last(), cells[i]);

        labelHashSet unique(indexList[i]);
        ASSERT_EQ(unique.size(), indexList[i].size());
    }
}
