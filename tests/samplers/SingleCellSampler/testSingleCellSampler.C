#include "fvCFD.H"
#include "SampledPGradField.H"
#include "SingleCellSampler.H"
#include "scalarListIOList.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "fixtures.H"


class SingleCellSamplerTest : public ChannelFlow
{};

TEST_F(SingleCellSamplerTest, ConstructorDefaults)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());
}

TEST_F(SingleCellSamplerTest, ConstructorCellCrawlingHIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "Crawling",
        "CubeRootVol",
        true
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), true);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
}


TEST_F(SingleCellSamplerTest, ConstructorCellCrawlingHIsDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "Crawling",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
    
}


TEST_F(SingleCellSamplerTest, ConstructorCellPointCrawlingHIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Crawling",
        "CubeRootVol",
        true
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), true);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.10);
    }
    
}

TEST_F(SingleCellSamplerTest, ConstructorCellPointCrawlingHIsDistance)
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

    h.boundaryFieldRef()[patch.index()] == 0.35;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Crawling",
        "WallNormalDistance",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.lengthScaleType(), "WallNormalDistance");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.35);
    }
    
}


TEST_F(SingleCellSamplerTest, ConstructorCellPointCrawlingHIsDistanceFirstCell)
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

    h.boundaryFieldRef()[patch.index()] == 0.05;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Crawling",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
    
}


TEST_F(SingleCellSamplerTest, ConstructorCellTreeHIsIndex)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    ASSERT_DEATH
    (
        {
            SingleCellSampler sampler
            (
                "SingleCellSampler",
                patch,
                3.0,
                "cell",
                "Tree",
                "CubeRootVol",
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerTest, ConstructorCellTreeHIsDistance)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "Tree",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
}


TEST_F(SingleCellSamplerTest, ConstructorCellPointTreeHIsDistance)
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

    h.boundaryFieldRef()[patch.index()] == 0.35;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Tree",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.35);
    }
}



TEST_F(SingleCellSamplerTest, ConstructorCellPointTreeHIsDistanceFirstCell)
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

    h.boundaryFieldRef()[patch.index()] == 0.05;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Tree",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
    
}


TEST_F(SingleCellSamplerTest, Sample)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    h.boundaryFieldRef()[patch.index()] == 3;

    auto & U = const_cast<volVectorField &>
    (
        mesh.thisDb().lookupObject<volVectorField>("U")
    );

    U.primitiveFieldRef() = vector(1, 0, 0);

    SingleCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        0.02,
        "cell",
        "Crawling",
        "CubeRootVol",
        3
    );
    
    U.primitiveFieldRef() = vector(4, 0, 0);
    
    sampler.sample();

    auto & sampledU = const_cast<scalarListIOList &>
    (
        sampler.db().lookupObject<scalarListIOList>("U")
    );
    
    forAll(sampledU, i)
    {
        ASSERT_FLOAT_EQ(sampledU[i][0], 2.5);
    }
}


TEST_F(SingleCellSamplerTest, AddField)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    h.boundaryFieldRef()[patch.index()] == 3;

    // auto sampledPGrad = SampledPGradField(patch);
    // auto sampledPGrad = SampledVelocityField(patch);

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        0.02,
        "cell",
        "Crawling",
        "CubeRootVol",
        3
    );
    
    sampler.addField(new SampledPGradField(patch));
    
    ASSERT_EQ(sampler.nSampledFields(), 3);
    ASSERT_TRUE(sampler.db().foundObject<scalarListIOList>("pGrad"));
}


TEST_F(SingleCellSamplerTest, createLengthListCubeRootVol)
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

    h.boundaryFieldRef()[patch.index()] == 0.5;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0
    );

    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    
    for (int i=0; i< sampler.lengthList().size(); i++)
    {
        ASSERT_FLOAT_EQ(sampler.lengthList()[i], 0.9283177667225558);
    }
}

TEST_F(SingleCellSamplerTest, createLengthListWallNormalDistance)
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

    h.boundaryFieldRef()[patch.index()] == 0.5;

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "Crawling",
        "WallNormalDistance"
    );

    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    
    for (int i=0; i< sampler.lengthList().size(); i++)
    {
        ASSERT_FLOAT_EQ(sampler.lengthList()[i], 0.2);
    }
}