#include "codeRules.H"
#include "fvCFD.H"
#include "SingleCellSampler.H"
#include "SampledVelocityField.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"


namespace Foam
{
    class DummySampler : public Sampler

    {
        public:
            TypeName("DummySampler");

            DummySampler
            (
                const fvPatch& p,
                scalar averagingTime,
                const word interpolationType,
                const word cellFinderType,
                const word lengthScaleType,
                bool hIsIndex,
                bool excludeWallAdjacent

            )
            :
                Sampler(p, averagingTime, interpolationType, cellFinderType,
                        lengthScaleType, hIsIndex, excludeWallAdjacent)
            {}

            DummySampler
            (
                const word & samplerName,
                const fvPatch & p,
                scalar averagingTime,
                const word interpolationType,
                const word cellFinderType,
                const word lengthScaleType,
                bool hIsIndex,
                bool excludeWallAdjacent
            )
            :
                Sampler(p, averagingTime, interpolationType, cellFinderType,
                        lengthScaleType, hIsIndex, excludeWallAdjacent)
            {}
            
            DummySampler(const DummySampler &) = default;

            void sample() const override
            {} 

            void createIndexList() override
            {}
            
        // Destructor
            virtual ~DummySampler()
            {}

    };

    defineTypeNameAndDebug(DummySampler, 0);
    addToRunTimeSelectionTable(Sampler, DummySampler, SamplerRTSTable);
}

class SamplerTest : public ChannelFlow
{};


TEST_F(SamplerTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "crawling",
        "CubeRootVol",
        false,
        false
    );

    ASSERT_EQ(&sampler.Sampler::patch(), &patch);
    ASSERT_EQ(sampler.Sampler::averagingTime(), 3.0);
    ASSERT_EQ(sampler.Sampler::interpolationType(), "cell");
    ASSERT_EQ(sampler.Sampler::cellFinderType(), "crawling");
    ASSERT_EQ(sampler.Sampler::lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler.Sampler::hIsIndex(), false);
    ASSERT_EQ(&sampler.Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler.Sampler::nSampledFields(), 0);
    ASSERT_EQ(sampler.Sampler::excludeWallAdjacent(), false);
    ASSERT_TRUE(mesh.foundObject<objectRegistry>("wallModelSampling"));
    ASSERT_TRUE
    (
        mesh.subRegistry("wallModelSampling").foundObject<objectRegistry>(patch.name())
    );
}


TEST_F(SamplerTest, NewNamePatchAveragingTime)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["bottomWall"];
    autoPtr<Sampler> sampler
    (
        Sampler::New
        (
            "DummySampler",
            patch,
            3.0,
            "cell",
            "crawling",
            "CubeRootVol",
            false,
            false
        )
    );

    ASSERT_EQ(sampler().type(), word("DummySampler"));
    ASSERT_EQ(&sampler().Sampler::patch(), &patch);
    ASSERT_EQ(sampler().Sampler::averagingTime(), 3.0);
    ASSERT_EQ(sampler().Sampler::interpolationType(), "cell");
    ASSERT_EQ(sampler().Sampler::cellFinderType(), "crawling");
    ASSERT_EQ(sampler().Sampler::lengthScaleType(), "CubeRootVol");
    ASSERT_EQ(sampler().Sampler::hIsIndex(), false);
    ASSERT_EQ(&sampler().Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler().Sampler::nSampledFields(), 0);
    ASSERT_EQ(sampler().Sampler::excludeWallAdjacent(), false);
    ASSERT_TRUE(mesh.foundObject<objectRegistry>("wallModelSampling"));
    ASSERT_TRUE
    (
        mesh.subRegistry("wallModelSampling").foundObject<objectRegistry>(patch.name())
    );
}


TEST_F(SamplerTest, NewDictionaryPatch)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    dictionary dict = dictionary();
    dict.lookupOrAddDefault(word("averagingTime"), 3.0);
    dict.lookupOrAddDefault(word("type"), word("DummySampler"));
    dict.lookupOrAddDefault(word("interpolationType"), word("cell"));
    dict.lookupOrAddDefault(word("sampler"), word("Crawling"));
    dict.lookupOrAddDefault(word("hIsIndex"), false);
    dict.lookupOrAddDefault(word("excludeWallAdjacent"), true);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    autoPtr<Sampler> sampler(Sampler::New(dict, patch));

    ASSERT_EQ(sampler().type(), word("DummySampler"));
    ASSERT_EQ(&sampler().Sampler::patch(), &patch);
    ASSERT_EQ(sampler().Sampler::averagingTime(), 3.0);
    ASSERT_EQ(sampler().Sampler::interpolationType(), "cell");
    ASSERT_EQ(sampler().Sampler::cellFinderType(), "Crawling");
    ASSERT_EQ(sampler().Sampler::hIsIndex(), false);
    ASSERT_EQ(&sampler().Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler().Sampler::nSampledFields(), 0);
    ASSERT_EQ(sampler().Sampler::excludeWallAdjacent(), true);
    ASSERT_TRUE(mesh.foundObject<objectRegistry>("wallModelSampling"));
    ASSERT_TRUE
    (
        mesh.subRegistry("wallModelSampling").foundObject<objectRegistry>(patch.name())
    );
}


TEST_F(SamplerTest, CreateFields)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "crawling",
        "CubeRootVol",
        false,
        false
    );

    ASSERT_TRUE(mesh.foundObject<volScalarField>("samplingCells"));
}


TEST_F(SamplerTest, Copy)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "crawling",
        "CubeRootVol",
        false,
        false
    );
    sampler.addField(new SampledVelocityField(patch));

    DummySampler sampler2(sampler);

    ASSERT_EQ(sampler.type(), sampler2.type());
    ASSERT_EQ(&sampler.patch(), &sampler2.patch());
    ASSERT_EQ(sampler.averagingTime(), sampler2.averagingTime());
    ASSERT_EQ(sampler.interpolationType(), sampler2.interpolationType());
    ASSERT_EQ(sampler.cellFinderType(), sampler2.cellFinderType());
    ASSERT_EQ(sampler.hIsIndex(), sampler2.hIsIndex());
    ASSERT_EQ(&sampler.mesh(), &sampler2.mesh());
    ASSERT_EQ(sampler.nSampledFields(), sampler2.nSampledFields());
    ASSERT_EQ(sampler.excludeWallAdjacent(), sampler2.excludeWallAdjacent());
}


//TEST_F(SamplerTest, Project)
//{
    //extern argList * mainArgs;
    //argList & args = *mainArgs;
    //Time runTime(Foam::Time::controlDictName, args);
    //autoPtr<fvMesh> meshPtr = createMesh(runTime);
    //const fvMesh & mesh = meshPtr();
    //createSamplingHeightField(mesh);

    //const fvPatch & patch = mesh.boundary()["bottomWall"];
    //DummySampler sampler("SingleCellSampler", patch, 3.0);

    //vectorField field(patch.size(), vector(1, 2, 3));
    //sampler.wrapProject(field);
    //forAll (field, i)
    //{
        //ASSERT_EQ(field[i], vector(1, 0, 3));
    //}
//}


TEST_F(SamplerTest, AddField)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "crawling",
        "CubeRootVol",
        false,
        false
    );
    
    sampler.addField(new SampledVelocityField(patch));
    ASSERT_EQ(sampler.Sampler::nSampledFields(), 1);
}

