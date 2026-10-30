#include "fvCFD.H"
#include "SampledPGradField.H"
#include "SingleCellSampler.H"
#include "scalarListIOList.H"
#include "gtest.h"
#undef Log
#include "fixtures.H"


class SingleCellSamplerParallelTest : public ChannelFlow
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


// Construct both samplers with the same settings and compare the results
static void compareToSingleCellSampler
(
    const fvMesh & mesh,
    const fvPatch & patch,
    const scalar h,
    const word interpolationType,
    const word cellFinderType,
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
    SingleCellSampler parallelSampler
    (
        patch,
        0,
        interpolationType,
        cellFinderType + "Parallel",
        lengthScaleType,
        hIsIndex
    );

    ASSERT_EQ(parallelSampler.interpolationType(), sampler.interpolationType());
    ASSERT_EQ(parallelSampler.indexList(), sampler.indexList());
    ASSERT_EQ(parallelSampler.sampledCells().size(), patch.size());
    ASSERT_EQ(parallelSampler.map().constructSize(), patch.size());

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(parallelSampler.h()[i], sampler.h()[i]);
        ASSERT_FLOAT_EQ
        (
            parallelSampler.lengthList()[i],
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
    parallelSampler.sample();
    const scalarListIOList & sampledU =
        parallelSampler.db().lookupObject<scalarListIOList>("U");

    forAll(sampledU, i)
    {
        forAll(sampledU[i], j)
        {
            ASSERT_FLOAT_EQ(sampledU[i][j], expected[i][j]);
        }
    }
}


TEST_F(SingleCellSamplerParallelTest, ConstructorTreeParallel)
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
        "TreeParallel"
    );

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(sampler.interpolationType(), "cell");
    ASSERT_EQ(sampler.cellFinderType(), "TreeParallel");
    ASSERT_EQ(sampler.lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.hIsIndex(), false);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());
    ASSERT_EQ(sampler.sampledCells().size(), patch.size());
    ASSERT_EQ(sampler.map().constructSize(), patch.size());
    ASSERT_EQ(sampler.nSampledFields(), 2);
    ASSERT_TRUE(sampler.db().foundObject<scalarListIOList>("U"));
    ASSERT_TRUE(sampler.db().foundObject<scalarListIOList>("wallGradU"));
}


TEST_F(SingleCellSamplerParallelTest, MatchesSingleCellSamplerCrawling)
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
                for (const scalar h : {0.0, 1.0, 3.0, 12.0})
                {
                    compareToSingleCellSampler
                    (
                        mesh,
                        patch,
                        h,
                        interpolationTypes[interpI],
                        "Crawling",
                        lengthScaleTypes[lengthI],
                        true
                    );
                }

                for (const scalar h : {0.0, 0.1, 0.45, 1.25, 5.0})
                {
                    compareToSingleCellSampler
                    (
                        mesh,
                        patch,
                        h,
                        interpolationTypes[interpI],
                        "Crawling",
                        lengthScaleTypes[lengthI],
                        false
                    );
                }
            }
        }
    }
}


TEST_F(SingleCellSamplerParallelTest, MatchesSingleCellSamplerTree)
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
                for (const scalar h : {0.0, 0.1, 0.45, 1.25, 5.0})
                {
                    compareToSingleCellSampler
                    (
                        mesh,
                        patch,
                        h,
                        interpolationTypes[interpI],
                        "Tree",
                        lengthScaleTypes[lengthI],
                        false
                    );
                }
            }
        }
    }
}


TEST_F(SingleCellSamplerParallelTest, HIsIndexFallsBackToCellInterpolation)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 2);

    SingleCellSampler sampler
    (
        patch,
        3.0,
        "cellPoint",
        "CrawlingParallel",
        "CubeRootVol",
        true
    );

    ASSERT_EQ(sampler.interpolationType(), "cell");

    forAll(sampler.h(), i)
    {
        ASSERT_FLOAT_EQ(sampler.h()[i], 0.3);
    }
}


TEST_F(SingleCellSamplerParallelTest, ConstructorTreeHIsIndex)
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
                patch,
                3.0,
                "cell",
                "TreeParallel",
                "CubeRootVol",
                true
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerParallelTest, AttemptInvalidLengthScaleName)
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
                patch,
                3.0,
                "cell",
                "TreeParallel",
                "RandomName"
            );
        },
        "FATAL ERROR"
    );
}


TEST_F(SingleCellSamplerParallelTest, SamplingCellsField)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 4);

    SingleCellSampler sampler
    (
        patch,
        3.0,
        "cell",
        "CrawlingParallel",
        "CubeRootVol",
        true
    );

    const volScalarField & samplingCells =
        mesh.lookupObject<volScalarField>("samplingCells");

    label nMarked = 0;

    forAll(samplingCells, cellI)
    {
        if (samplingCells[cellI] == patch.index())
        {
            nMarked++;
            ASSERT_FLOAT_EQ(mesh.C()[cellI][1], 0.7);
        }
    }

    ASSERT_EQ(nMarked, patch.size());
}


TEST_F(SingleCellSamplerParallelTest, Sample)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 3);

    auto & U = const_cast<volVectorField &>
    (
        mesh.thisDb().lookupObject<volVectorField>("U")
    );

    U.primitiveFieldRef() = vector(1, 0, 0);

    SingleCellSampler sampler
    (
        patch,
        0.02,
        "cell",
        "CrawlingParallel",
        "CubeRootVol",
        true
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


TEST_F(SingleCellSamplerParallelTest, AddField)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);
    createPGradField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 3);

    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");
    pGrad.primitiveFieldRef() = mesh.C();

    SingleCellSampler sampler
    (
        patch,
        0,
        "cell",
        "CrawlingParallel",
        "CubeRootVol",
        true
    );

    sampler.addField(new SampledPGradField(patch));

    ASSERT_EQ(sampler.nSampledFields(), 3);
    ASSERT_TRUE(sampler.db().foundObject<scalarListIOList>("pGrad"));

    // Recomputing would overwrite pGrad with the gradient of p
    sampler.sample();

    const scalarListIOList & sampledPGrad =
        sampler.db().lookupObject<scalarListIOList>("pGrad");

    forAll(sampledPGrad, i)
    {
        ASSERT_FLOAT_EQ(sampledPGrad[i][0], pGrad[sampler[i]][0]);
        ASSERT_EQ(sampledPGrad[i][1], 0);
        ASSERT_FLOAT_EQ(sampledPGrad[i][2], pGrad[sampler[i]][2]);
    }
}


TEST_F(SingleCellSamplerParallelTest, Copy)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 3);

    SingleCellSampler sampler
    (
        patch,
        3.0,
        "cell",
        "CrawlingParallel",
        "WallNormalDistance",
        true
    );

    SingleCellSampler copied(sampler);

    ASSERT_EQ(copied.cellFinderType(), "CrawlingParallel");
    ASSERT_EQ(copied.indexList(), sampler.indexList());
    ASSERT_EQ(copied.sampledCells(), sampler.sampledCells());
    ASSERT_EQ(copied.h(), sampler.h());
    ASSERT_EQ(copied.lengthList(), sampler.lengthList());
    ASSERT_EQ(copied.nSampledFields(), sampler.nSampledFields());
    ASSERT_EQ(copied.map().constructSize(), sampler.map().constructSize());
    ASSERT_EQ(copied.map().subMap(), sampler.map().subMap());
    ASSERT_EQ(copied.map().constructMap(), sampler.map().constructMap());
}


TEST_F(SingleCellSamplerParallelTest, Write)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    setH(mesh, patch, 3);

    autoPtr<SingleCellSampler> sampler
    (
        new SingleCellSampler
        (
            patch,
            3.0,
            "cell",
            "CrawlingParallel",
            "CubeRootVol",
            true
        )
    );

    OStringStream os;
    sampler->write(os);

    const std::string output(os.str().c_str());

    ASSERT_NE(output.find("CrawlingParallel;"), std::string::npos);
    ASSERT_NE(output.find("interpolationType"), std::string::npos);
    ASSERT_NE(output.find("hIsIndex"), std::string::npos);
}
