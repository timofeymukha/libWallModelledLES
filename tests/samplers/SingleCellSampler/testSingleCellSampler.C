#include "fvCFD.H"
#include "SingleCellSampler.H"
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
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());
}

TEST_F(SingleCellSamplerTest, ConstructorCellCrawlingTrue)
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
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
}


TEST_F(SingleCellSamplerTest, ConstructorCellCrawlingFalse)
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
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
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


TEST_F(SingleCellSamplerTest, ConstructorCellPointCrawlingTrue)
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
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.10);
    }
    
}

TEST_F(SingleCellSamplerTest, ConstructorCellPointCrawlingFalse)
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
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cellPoint",
        "Crawling",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Crawling");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.11);
    }
    
}


TEST_F(SingleCellSamplerTest, ConstructorCellTreeTrue)
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
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerTest, ConstructorCellTreeFalse)
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
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.1);
    }
    
    // const scalarField & hPatch = h.boundaryFieldRef()[patch.index()] == 10.0;

}


TEST_F(SingleCellSamplerTest, ConstructorCellPointTreeFalse)
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
        "Tree",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "Tree");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.11);
    }
}