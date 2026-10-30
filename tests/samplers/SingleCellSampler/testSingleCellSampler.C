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

// Set the value of h at a given patch
static void setH(const fvMesh & mesh, const fvPatch & patch, const scalar value)
{
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    h.boundaryFieldRef()[patch.index()] == value;
}


// Construct samplers with the same settings, except for the cell finder,
// and compare the results
static void compareCellFinders
(
    const fvMesh & mesh,
    const fvPatch & patch,
    const scalar h,
    const word interpolationType,
    const word cellFinderType,
    const word otherCellFinderType,
    const word lengthScaleType,
    const bool hIsIndex
)
{
    setH(mesh, patch, h);
    SingleCellSampler sampler
    (
        patch,
        0,
        interpolationType,
        cellFinderType,
        lengthScaleType,
        hIsIndex
    );

    setH(mesh, patch, h);
    SingleCellSampler otherSampler
    (
        patch,
        0,
        interpolationType,
        otherCellFinderType,
        lengthScaleType,
        hIsIndex
    );

    ASSERT_EQ(otherSampler.interpolationType(), sampler.interpolationType());
    ASSERT_EQ(otherSampler.indexList(), sampler.indexList());
    ASSERT_EQ(otherSampler.sampledCells().size(), patch.size());
    ASSERT_EQ(otherSampler.map().constructSize(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(otherSampler.h()[i], sampler.h()[i]);
        ASSERT_FLOAT_EQ
        (
            otherSampler.lengthList()[i],
            sampler.lengthList()[i]
        );
    }

    // The sampled values, both samplers store them in the same registry
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");
    U.primitiveFieldRef() = mesh.C();
    U.boundaryFieldRef()[patch.index()] == vector(0.1, 0, 0.2);

    sampler.sample();
    const scalarListList expected
    (
        sampler.db().lookupObject<scalarListIOList>("U")
    );

    U.primitiveFieldRef() = vector::zero;
    sampler.sample();

    U.primitiveFieldRef() = mesh.C();
    otherSampler.sample();
    const scalarListIOList & sampledU =
        otherSampler.db().lookupObject<scalarListIOList>("U");

    forAll(sampledU, i)
    {
        forAll(sampledU[i], j)
        {
            ASSERT_FLOAT_EQ(sampledU[i][j], expected[i][j]);
        }
    }
}


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
    ASSERT_EQ(sampler.cellFinderType(), "NormalRay");
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


TEST_F(SingleCellSamplerTest, ConstructorCellNormalRayHIsIndex)
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
                "NormalRay",
                "CubeRootVol",
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerTest, AttemptInvalidLengthScaleName)
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
                "NormalRay",
                "RandomName",
                false
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerTest, ConstructorCellNormalRayHIsDistance)
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
        "NormalRay",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "NormalRay");
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


TEST_F(SingleCellSamplerTest, ConstructorCellPointNormalRayHIsDistance)
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
        "NormalRay",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "NormalRay");
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



TEST_F(SingleCellSamplerTest, ConstructorCellPointNormalRayHIsDistanceFirstCell)
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
        "NormalRay",
        "CubeRootVol",
        false
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cellPoint");
    ASSERT_EQ(sampler.cellFinderType(), "NormalRay");
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


TEST_F(SingleCellSamplerTest, SampledCellsAndMap)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 0.45);

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "NormalRay"
    );

    // In serial, all sampling cells are on this processor
    ASSERT_EQ(sampler.sampledCells().size(), patch.size());
    ASSERT_EQ(sampler.map().constructSize(), patch.size());

    labelList cells(sampler.sampledCells());
    sampler.map().distribute(cells);
    ASSERT_EQ(cells, sampler.indexList());

    forAll(cells, i)
    {
        ASSERT_FLOAT_EQ(mesh.C()[cells[i]][1], 0.5);
    }
}


// On the layered channel mesh, the Crawling and NormalRay cell finders select
// the same cells for distances
TEST_F(SingleCellSamplerTest, CrawlingMatchesNormalRay)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    const wordList patches({"bottomWall", "topWall"});
    const wordList interpolationTypes({"cell", "cellPoint"});
    const wordList lengthScaleTypes({"CubeRootVol", "WallNormalDistance"});

    forAll(patches, patchI)
    {
        const fvPatch & patch = mesh.boundary()[patches[patchI]];

        forAll(interpolationTypes, interpI)
        {
            forAll(lengthScaleTypes, lengthI)
            {
                for (const scalar h : {0.0, 0.1, 0.45, 1.25, 1.9})
                {
                    compareCellFinders
                    (
                        mesh,
                        patch,
                        h,
                        interpolationTypes[interpI],
                        "Crawling",
                        "NormalRay",
                        lengthScaleTypes[lengthI],
                        false
                    );
                }
            }
        }
    }
}
