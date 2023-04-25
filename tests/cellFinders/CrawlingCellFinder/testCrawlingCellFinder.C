#include "fvCFD.H"
#include "CrawlingCellFinder.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class CrawlingCellFinderTest : public ChannelFlow
{};

TEST_F(CrawlingCellFinderTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const fvPatch & patch = mesh.boundary()["bottomWall"];
    h.boundaryFieldRef()[patch.index()] == 2;
    CrawlingCellFinder finder(patch);
}

TEST_F(CrawlingCellFinderTest, FindIndexBasedBottomWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }

    h.boundaryFieldRef()[patch.index()] == 2;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.3);
    }

    h.boundaryFieldRef()[patch.index()] == 10;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }
}

TEST_F(CrawlingCellFinderTest, FindIndexBasedTopWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }

    h.boundaryFieldRef()[patch.index()] == 2;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.7);
    }

    h.boundaryFieldRef()[patch.index()] == 10;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, FindIndexBasedInvalidIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == -1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }
}


TEST_F(CrawlingCellFinderTest, FindIndexBasedTooLargeIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 100;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceBasedBottomWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }

    h.boundaryFieldRef()[patch.index()] == 0.32;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.3);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceBasedTopWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }

    h.boundaryFieldRef()[patch.index()] == 0.32;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.7);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceInvalidDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == -0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}

TEST_F(CrawlingCellFinderTest, FindDistanceZeroDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.0;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceTooLargeDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 10.0;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }
}


TEST_F(CrawlingCellFinderTest, FindIndexBasedBottomWallMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    h.boundaryFieldRef()[patch.index()] == 2;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 2);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
    }
}


TEST_F(CrawlingCellFinderTest, FindIndexBasedTopWallMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 1.9);
    }

    h.boundaryFieldRef()[patch.index()] == 5;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 5);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 1.9);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 1.7);
        ASSERT_FLOAT_EQ(C[indexList[i][4]][1], 1.1);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceBasedBottomWallMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    h.boundaryFieldRef()[patch.index()] == 0.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 5);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
        ASSERT_FLOAT_EQ(C[indexList[i][4]][1], 0.9);
    }
}


TEST_F(CrawlingCellFinderTest, FindDistanceBasedTopWallMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["topWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 1.9);
    }

    h.boundaryFieldRef()[patch.index()] == 0.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 5);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 1.9);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 1.7);
        ASSERT_FLOAT_EQ(C[indexList[i][4]][1], 1.1);
    }
}


TEST_F(CrawlingCellFinderTest, NegativeIndexMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == -5;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, NegativeDistanceMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == -5;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }
}


TEST_F(CrawlingCellFinderTest, LargeIndexMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 555;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 10);
        ASSERT_FLOAT_EQ(C[indexList[i][9]][1], 1.9);
    }
}


TEST_F(CrawlingCellFinderTest, LargeDistanceMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 555;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 10);
        ASSERT_FLOAT_EQ(C[indexList[i][9]][1], 1.9);
    }
}


TEST_F(CrawlingCellFinderTest, ExcludeWallAdjacentMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 4;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        true
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 3);
        ASSERT_FLOAT_EQ(C[indexList[i][2]][1], 0.7);
    }
}


TEST_F(CrawlingCellFinderTest, ExcludeWallAdjacentOneCellMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    const volVectorField & C = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    CrawlingCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true,
        true
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }
}