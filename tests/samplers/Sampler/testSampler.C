#include "codeRules.H"
#include "fvCFD.H"
#include "SingleCellSampler.H"
#include "SampledVelocityField.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"


namespace Foam
{
    class DummySampler : public Sampler

    {
        public:
            TypeName("DummySampler");

            DummySampler
            (
                const fvPatch& p,
                scalar averagingTime
            )
            :
                Sampler(p, averagingTime)
            {}

            DummySampler
            (
                const word & samplerName,
                const fvPatch & p,
                scalar averagingTime
            )
            :
                Sampler(samplerName, p, averagingTime)
            {}
            
            DummySampler(const DummySampler &) = default;

            void sample() const override
            {} 

            void createIndexList() override
            {}
            
        // Destructor
            virtual ~DummySampler()
            {}

            void wrapProject(vectorField & field) const
            {
                return Sampler::project(field);
            }

            tmp<volScalarField> wrapDistanceField() const
            {
                return Sampler::distanceField();
            }

            tmp<labelField> wrapFindSearchCellLabels() const
            {
                return Sampler::findSearchCellLabels();
            }
            
    };

    defineTypeNameAndDebug(DummySampler, 0);
    addToRunTimeSelectionTable(Sampler, DummySampler, PatchAndAveragingTime);
}

class SamplerTest : public ::testing::Test
{


    public:
        SamplerTest()
        {
            system("cp -r testCases/channel_flow/system .");
            system("cp -r testCases/channel_flow/constant .");
            system("cp -r testCases/channel_flow/0 .");

        }

        ~SamplerTest()
        {
            system("rm -r system");
            system("rm -r constant");
            system("rm -r 0");
        }
};


void createSamplingHeightField(const fvMesh & mesh)
{
    mesh.time().store
    (
        new volScalarField
        (
            IOobject
            (
                "h",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}


TEST_F(SamplerTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    ASSERT_EQ(&sampler.Sampler::patch(), &patch);
    ASSERT_EQ(sampler.Sampler::averagingTime(), 3.0);
    ASSERT_EQ(&sampler.Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler.Sampler::nSampledFields(), 0);
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
#include "createMesh.H"
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["bottomWall"];
    autoPtr<Sampler> sampler(Sampler::New("DummySampler", patch, 3.0));

    ASSERT_EQ(sampler().type(), word("DummySampler"));
    ASSERT_EQ(&sampler().Sampler::patch(), &patch);
    ASSERT_EQ(sampler().Sampler::averagingTime(), 3.0);
    ASSERT_EQ(&sampler().Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler().Sampler::nSampledFields(), 0);
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
#include "createMesh.H"
    createSamplingHeightField(mesh);

    dictionary dict = dictionary();
    dict.lookupOrAddDefault(word("averagingTime"), 3.0);
    dict.lookupOrAddDefault(word("type"), word("DummySampler"));

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    autoPtr<Sampler> sampler(Sampler::New(dict, patch));

    ASSERT_EQ(sampler().type(), word("DummySampler"));
    ASSERT_EQ(&sampler().Sampler::patch(), &patch);
    ASSERT_EQ(sampler().Sampler::averagingTime(), 3.0);
    ASSERT_EQ(&sampler().Sampler::mesh(), &mesh);
    ASSERT_EQ(sampler().Sampler::nSampledFields(), 0);
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
#include "createMesh.H"
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    ASSERT_TRUE(mesh.foundObject<volScalarField>("samplingCells"));
}


TEST_F(SamplerTest, Copy)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);
    sampler.addField(new SampledVelocityField(patch));

    DummySampler sampler2(sampler);

    ASSERT_EQ(sampler.type(), sampler2.type());
    ASSERT_EQ(&sampler.patch(), &sampler2.patch());
    ASSERT_EQ(sampler.averagingTime(), sampler2.averagingTime());
    ASSERT_EQ(&sampler.mesh(), &sampler2.mesh());
    ASSERT_EQ(sampler.nSampledFields(), sampler2.nSampledFields());
}


TEST_F(SamplerTest, Project)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    vectorField field(patch.size(), vector(1, 2, 3));
    sampler.wrapProject(field);
    forAll (field, i)
    {
        ASSERT_EQ(field[i], vector(1, 0, 3));
    }
}


TEST_F(SamplerTest, AddField)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);
    
    sampler.addField(new SampledVelocityField(patch));
    ASSERT_EQ(sampler.Sampler::nSampledFields(), 1);
}


TEST_F(SamplerTest, DistanceFieldBottom)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    tmp<volScalarField> tDist = sampler.wrapDistanceField();
    const volScalarField & dist = tDist();

    const vectorField C = mesh.C();
    
    forAll (C, i)
    {
        ASSERT_NEAR(dist[i], C[i][1], 1e-8);
    }
}

TEST_F(SamplerTest, DistanceFieldTop)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    createSamplingHeightField(mesh);


    const fvPatch & patch = mesh.boundary()["topWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    tmp<volScalarField> tDist = sampler.wrapDistanceField();
    const volScalarField & dist = tDist();

    const vectorField C = mesh.C();
    
    forAll (C, i)
    {
        ASSERT_NEAR(dist[i], 2 - C[i][1], 1e-8);
    }
}


TEST_F(SamplerTest, FindSearchCellLabels)
{
    extern argList * mainArgs;
    argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"
    //system("funkySetFields -field h -expression \"0.0\" -latestTime -valuePatches \"bottomWall topWall\"");
    createSamplingHeightField(mesh);
    volScalarField & h = 
        const_cast<volScalarField &>(mesh.lookupObject<volScalarField> ("h"));
    const fvPatch & patchBottom = mesh.boundary()["topWall"];

    h.boundaryFieldRef()[patchBottom.index()] == 0;
    DummySampler samplerBottom("SingleCellSampler", patchBottom, 3.0);
    tmp<labelField> tCellLabelsBottom = samplerBottom.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsBottom().size(), 0);

    h.boundaryFieldRef()[patchBottom.index()] == 0.1;
    tCellLabelsBottom = samplerBottom.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsBottom().size(), 9);

    h.boundaryFieldRef()[patchBottom.index()] == 0.21;
    tCellLabelsBottom = samplerBottom.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsBottom().size(), 18);

    h.boundaryFieldRef()[patchBottom.index()] == 1;
    tCellLabelsBottom = samplerBottom.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsBottom().size(), 90);

    const fvPatch & patchTop = mesh.boundary()["topWall"];

    h.boundaryFieldRef()[patchTop.index()] == 0;
    DummySampler samplerTop("SingleCellSampler", patchTop, 3.0);
    tmp<labelField> tCellLabelsTop = samplerTop.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsTop().size(), 0) << "h = 0";

    h.boundaryFieldRef()[patchTop.index()] == 0.1;
    tCellLabelsTop = samplerTop.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsTop().size(), 9) << "h = 0.1";

    h.boundaryFieldRef()[patchTop.index()] == 0.21;
    tCellLabelsTop = samplerTop.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsTop().size(), 18) << "h = 0.21";

    h.boundaryFieldRef()[patchTop.index()] == 1;
    tCellLabelsTop = samplerTop.wrapFindSearchCellLabels();
    ASSERT_EQ(tCellLabelsTop().size(), 90) << "h = 1";
    
}
