#include "fvCFD.H"
#include "CrawlingCellFinder.H"
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        mesh.thisDb().lookupObject<volScalarField>("h")
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


TEST_F(CrawlingCellFinderTest, MulticellFindIndexBasedBottomWall)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("h")
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
        true
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
        true
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 2);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
    }
}