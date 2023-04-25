#include "fvCFD.H"
#include "TreeCellFinder.H"
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"
#include "fixtures.H"


class TreeCellFinderTest : public ChannelFlow
{};

TEST_F(TreeCellFinderTest, FullConstructor)
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
    h.boundaryFieldRef()[patch.index()] == 0.25;
    TreeCellFinder finder(patch);
}


TEST_F(TreeCellFinderTest, FindDistanceBasedBottomWallSingleCell)
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
    TreeCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }

    h.boundaryFieldRef()[patch.index()] == 0.32;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.3);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }
}


TEST_F(TreeCellFinderTest, FindDistanceBasedTopWallSingleCell)
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
    TreeCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.9);
    }

    h.boundaryFieldRef()[patch.index()] == 0.32;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 1.7);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(TreeCellFinderTest, FindDistanceInvalidDistanceSingleCell)
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
    TreeCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == -0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}

TEST_F(TreeCellFinderTest, FindDistanceZeroDistanceSingleCell)
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
    TreeCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.0;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(TreeCellFinderTest, FindDistanceTooLargeDistanceSingleCell)
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
    TreeCellFinder finder(patch);
    labelList indexList(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 10.0;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexList, i)
    {
        ASSERT_FLOAT_EQ(C[indexList[i]][1], 0.1);
    }
}


TEST_F(TreeCellFinderTest, DistanceFieldBottom)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["bottomWall"];
    TreeCellFinder finder(patch);

    tmp<volScalarField> tDist = finder.distanceField();
    const volScalarField & dist = tDist();

    const vectorField C = mesh.C();
    
    forAll (C, i)
    {
        ASSERT_NEAR(dist[i], C[i][1], 1e-8);
    }
}

TEST_F(TreeCellFinderTest, DistanceFieldTop)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["topWall"];
    TreeCellFinder finder(patch);

    tmp<volScalarField> tDist = finder.distanceField();
    const volScalarField & dist = tDist();

    const vectorField C = mesh.C();
    
    forAll (C, i)
    {
        ASSERT_NEAR(dist[i], 2 - C[i][1], 1e-8);
    }
}


TEST_F(TreeCellFinderTest, FindCandidateCellLabelsBottom)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patchBottom = mesh.boundary()["topWall"];
    TreeCellFinder finderBottom(patchBottom);

    scalarField h(finderBottom.patch().size());

    h = 0.0;
    tmp<Foam::volScalarField> distField(finderBottom.distanceField());

    const scalarField & dist =
#ifdef FOAM_NEW_GEOMFIELD_RULES
        distField().primitiveField();
#else
        distField().internalField();
#endif
    tmp<labelField> tCellLabelsBottom =
        finderBottom.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsBottom().size(), 0);

    h = 0.1;
    tCellLabelsBottom = finderBottom.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsBottom().size(), 9);

    h = 0.21;
    tCellLabelsBottom = finderBottom.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsBottom().size(), 18);

    h = 1;
    tCellLabelsBottom = finderBottom.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsBottom().size(), 90);

}


TEST_F(TreeCellFinderTest, FindCandidateCellLabelsTop)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patchTop = mesh.boundary()["topWall"];
    TreeCellFinder finderTop(patchTop);

    scalarField h(finderTop.patch().size());

    tmp<Foam::volScalarField> distField = finderTop.distanceField();
    const scalarField & dist =
#ifdef FOAM_NEW_GEOMFIELD_RULES
        distField().primitiveField();
#else
        distField().internalField();
#endif

    h = 0.0;
    tmp<labelField> tCellLabelsTop = finderTop.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsTop().size(), 0) << "h = 0";

    h = 0.1;
    tCellLabelsTop = finderTop.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsTop().size(), 9) << "h = 0.1";

    h = 0.21;
    tCellLabelsTop = finderTop.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsTop().size(), 18) << "h = 0.21";

    h = 1;
    tCellLabelsTop = finderTop.findCandidateCellLabels(dist, h);
    ASSERT_EQ(tCellLabelsTop().size(), 90) << "h = 1";
}

TEST_F(TreeCellFinderTest, FindDistanceBasedBottomWallMultiCell)
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
    TreeCellFinder finder(patch);
    labelListList indexList(patch.size());

    // Single, wall-adjacent cell
    h.boundaryFieldRef()[patch.index()] == 0.1;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    // h is boundary between first two cells, expecting a 1-cell list
    h.boundaryFieldRef()[patch.index()] == 0.2;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 1);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
    }

    // h above 3rd cell center, expecting a 3-cell list
    h.boundaryFieldRef()[patch.index()] == 0.55;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 3);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.3);
        ASSERT_FLOAT_EQ(C[indexList[i][2]][1], 0.5);
    }
}

TEST_F(TreeCellFinderTest, FindDistanceBasedBottomWallWithExclusionMultiCell)
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
    TreeCellFinder finder(patch);
    labelListList indexList(patch.size());

    // Single, wall-adjacent cell
    h.boundaryFieldRef()[patch.index()] == 0.1;
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

    // h is boundary between first two cells, expecting a 1-cell list
    h.boundaryFieldRef()[patch.index()] == 0.2;
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

    // h above 3rd cell center, expecting a 2-cell list
    h.boundaryFieldRef()[patch.index()] == 0.55;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 2);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.3);
        ASSERT_FLOAT_EQ(C[indexList[i][1]][1], 0.5);
    }
}

TEST_F(TreeCellFinderTest, ZeroDistanceMultiCell)
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
    TreeCellFinder finder(patch);
    labelListList indexList(patch.size());

    // Expect wall-adjacent
    h.boundaryFieldRef()[patch.index()] == 0.0;
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
}

TEST_F(TreeCellFinderTest, InvalidDistanceMultiCell)
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
    TreeCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == -10.0;
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
}

TEST_F(TreeCellFinderTest, TooLargeDistanceMultiCell)
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
    TreeCellFinder finder(patch);
    labelListList indexList(patch.size());

    h.boundaryFieldRef()[patch.index()] == 2.9;
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexList, i)
    {
        ASSERT_EQ(indexList[i].size(), 10);
        ASSERT_FLOAT_EQ(C[indexList[i][0]][1], 0.1);
        ASSERT_FLOAT_EQ(C[indexList[i][9]][1], 1.9);
    }
}