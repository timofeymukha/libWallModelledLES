#include "fvCFD.H"
#include "CrawlingCellFinder.H"
#include "NormalRayCellFinder.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class CellFinderCompatibilityTest : public ChannelFlow
{};


TEST_F(CellFinderCompatibilityTest, ExactDistance)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListNormalRay(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 0.3;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

}


TEST_F(CellFinderCompatibilityTest, NegativeDistance)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListNormalRay(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == -0.1;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

}


// Both use the last cell before the opposite wall
TEST_F(CellFinderCompatibilityTest, TooLargeDistance)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListNormalRay(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 3.1;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

}


TEST_F(CellFinderCompatibilityTest, InexactDistance)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListNormalRay(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.05;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 0.15;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.85;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.95;
    findCellIndices
    (
        finderCrawling,
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    findCellIndices
    (
        finderNormalRay,
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListNormalRay[i]);
    }
}


TEST_F(CellFinderCompatibilityTest, ExactDistanceMulticell)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListNormalRay(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

    h.boundaryFieldRef()[patch.index()] == 0.3;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }
}

TEST_F(CellFinderCompatibilityTest, ExactDistanceExcludeMulticell)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListNormalRay(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

    h.boundaryFieldRef()[patch.index()] == 0.3;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true 
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }
}


TEST_F(CellFinderCompatibilityTest, NegativeDistanceMulticell)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListNormalRay(patch.size());

    h.boundaryFieldRef()[patch.index()] == -1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

}


TEST_F(CellFinderCompatibilityTest, TooLargeDistanceMulticell)
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

    CrawlingCellFinder finderCrawling(patch);
    NormalRayCellFinder finderNormalRay(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListNormalRay(patch.size());

    h.boundaryFieldRef()[patch.index()] == 10;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderNormalRay.findCellIndices
    (
        indexListNormalRay,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListNormalRay[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListNormalRay[i][j]);
        }
    }

}
