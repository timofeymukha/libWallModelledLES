#include "fvCFD.H"
#include "SingleCellSampler.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"


class SamplerTest : public ::testing::Test
{
    public:
        SamplerTest()
        {
            system("cp -r ../../testCases/channel_flow/system .");
            system("cp -r ../../testCases/channel_flow/constant .");
            system("cp -r ../../testCases/channel_flow/0 .");
        }

        ~SamplerTest()
        {
            system("rm -r system");
            system("rm -r constant");
            system("rm -r 0");
        }
};

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
            
    };

    defineTypeNameAndDebug(DummySampler, 0);
    addToRunTimeSelectionTable(Sampler, DummySampler, PatchAndAveragingTime);
}

TEST_F(SamplerTest, FullConstructor)
{
    extern int mainArgc;
    extern char** mainArgv;

    argList * argsPtr = new argList(mainArgc, mainArgv);
    argList & args = *argsPtr;
    Time runTime(Foam::Time::controlDictName, args);
#include "createMesh.H"

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

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    DummySampler sampler("SingleCellSampler", patch, 3.0);

    ASSERT_EQ(&sampler.Sampler::patch(), &patch);
    ASSERT_EQ(sampler.Sampler::averagingTime(), 3.0);
    ASSERT_EQ(&sampler.Sampler::mesh(), &mesh);

    EXPECT_EXIT(delete argsPtr, ::testing::ExitedWithCode(0), "");
}

int mainArgc;
char** mainArgv;

int main(int argc, char **argv)
{
    Foam::argList::noBanner();

    ::testing::InitGoogleMock(&argc, argv);
    ::testing::InitGoogleTest(&argc, argv);

    mainArgc = argc;
    mainArgv = argv;

    return RUN_ALL_TESTS();
}
