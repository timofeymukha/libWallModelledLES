#include "fvCFD.H"
#include "CrawlingCellFinder.H"
#include "TreeCellFinder.H"
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
    TreeCellFinder finderTree(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListTree(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 0.3;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.9;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
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
    TreeCellFinder finderTree(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListTree(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == -0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

}


// WE EXPECT DIFFERENCE BEHAVIOUR IN THIS CASE
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
    TreeCellFinder finderTree(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListTree(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 3.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_NE(indexListCrawling[i], indexListTree[i]);
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
    TreeCellFinder finderTree(patch);

    labelList indexListCrawling(patch.size(), 0);
    labelList indexListTree(patch.size(), 0);

    h.boundaryFieldRef()[patch.index()] == 0.05;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 0.15;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.85;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
    }

    h.boundaryFieldRef()[patch.index()] == 1.95;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()]
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_FLOAT_EQ(indexListCrawling[i], indexListTree[i]);
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
    TreeCellFinder finderTree(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListTree(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        false
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        false
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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
    TreeCellFinder finderTree(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListTree(patch.size());

    h.boundaryFieldRef()[patch.index()] == 0.1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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
    TreeCellFinder finderTree(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListTree(patch.size());

    h.boundaryFieldRef()[patch.index()] == -1;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
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
    TreeCellFinder finderTree(patch);

    labelListList indexListCrawling(patch.size());
    labelListList indexListTree(patch.size());

    h.boundaryFieldRef()[patch.index()] == 10;
    finderCrawling.findCellIndices
    (
        indexListCrawling,
        h.boundaryFieldRef()[patch.index()],
        false,
        true
    );

    finderTree.findCellIndices
    (
        indexListTree,
        h.boundaryFieldRef()[patch.index()],
        true
    );

    forAll(indexListCrawling, i)
    {
        ASSERT_EQ(indexListCrawling[i].size(), indexListTree[i].size());
        forAll(indexListCrawling[i], j)
        {
            ASSERT_FLOAT_EQ(indexListCrawling[i][j], indexListTree[i][j]);
        }
    }

}