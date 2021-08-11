#include "fvCFD.H"
#include "MultiCellSampler.H"
#include "SampledPGradField.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "fixtures.H"


class MultiCellSamplerTest : public ChannelFlow
{};

TEST_F(MultiCellSamplerTest, ConstructorDefaults)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    MultiCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        3.0
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());
}

TEST_F(MultiCellSamplerTest, ConstructorCellCrawling)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("h")
    );

    h.boundaryFieldRef()[patch.index()] == 2;
    MultiCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        3.0,
        "cell",
        "Crawling",
        true
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.hIsIndex(), true);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_EQ(sampler.h()[i].size(), 2);
        ASSERT_FLOAT_EQ(sampler.h()[i][1], 0.3);
    }
}


TEST_F(MultiCellSamplerTest, ConstructorAttemptCellPoint)
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
            MultiCellSampler sampler
            (
                "MultiCellSampler",
                patch,
                3.0,
                "cellPoint",
                "Crawling",
                false
            );
        },
        "FATAL ERROR"
    );
    
}


TEST_F(MultiCellSamplerTest, ConstructorAttemptTreeHIsIndex)
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
            MultiCellSampler sampler
            (
                "MultiCellSampler",
                patch,
                3.0,
                "cell",
                "Tree",
                true
            );
        },
        "FATAL ERROR"
    );

    
}


TEST_F(MultiCellSamplerTest, ConstructorTree)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("h")
    );

    h.boundaryFieldRef()[patch.index()] == 0.5;
    MultiCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        3.0,
        "cell",
        "Tree",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_EQ(sampler.h()[i].size(), 3);
        ASSERT_FLOAT_EQ(sampler.h()[i][0], 0.1);
        ASSERT_FLOAT_EQ(sampler.h()[i][1], 0.3);
        ASSERT_FLOAT_EQ(sampler.h()[i][2], 0.5);
    }
}


TEST_F(MultiCellSamplerTest, AttemptInvalidCellFinderName)
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
            MultiCellSampler sampler
            (
                "MultiCellSampler",
                patch,
                3.0,
                "cell",
                "RandomName",
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(MultiCellSamplerTest, Sample)
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
        mesh.thisDb().lookupObject<volScalarField>("h")
    );

    h.boundaryFieldRef()[patch.index()] == 3;

    auto & U = const_cast<volVectorField &>
    (
        mesh.thisDb().lookupObject<volVectorField>("U")
    );

    U.primitiveFieldRef() = vector(1, 0, 0);

    MultiCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        0.02,
        "cell",
        "Crawling",
        3
    );
    
    U.primitiveFieldRef() = vector(4, 0, 0);
    
    sampler.sample();

    auto & sampledU = const_cast<scalarListListIOList &>
    (
        sampler.db().lookupObject<scalarListListIOList>("U")
    );
    
    forAll(sampledU, i)
    {
        forAll(sampledU[i], j)
        {
            ASSERT_FLOAT_EQ(sampledU[i][j][0], 2.5);
        }
    }
}


TEST_F(MultiCellSamplerTest, AddField)
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
        mesh.thisDb().lookupObject<volScalarField>("h")
    );

    h.boundaryFieldRef()[patch.index()] == 3;

    MultiCellSampler sampler
    (
        "MultiCellSampler",
        patch,
        0.02,
        "cell",
        "Crawling",
        3
    );
    
    sampler.addField(new SampledPGradField(patch));
    
    ASSERT_EQ(sampler.nSampledFields(), 3);
    ASSERT_TRUE(sampler.db().foundObject<scalarListListIOList>("pGrad"));
}