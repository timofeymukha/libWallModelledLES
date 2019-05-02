#include "fvCFD.H"
#include "SingleCellSampler.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"


class SingleCellSamplerTest : public ::testing::Test
{
    public:

        //argList * argsPtr;
        //autoPtr<Time> runTimePtr;
        //autoPtr<fvMesh> meshPtr;

        SingleCellSamplerTest()
        //:
            //argsPtr(nullptr),
            //runTimePtr(nullptr),
            //meshPtr(nullptr)
        {
            //extern int mainArgc;
            //extern char** mainArgv;

            system("cp -r ../../testCases/channel_flow/system .");
            system("cp -r ../../testCases/channel_flow/constant .");
            system("cp -r ../../testCases/channel_flow/0 .");

            //argsPtr = new argList(mainArgc, mainArgv);
            //runTimePtr.reset(new Time(Foam::Time::controlDictName, *argsPtr));
            //Time & runTime = runTimePtr();

            //meshPtr.reset
            //(
                //new fvMesh
                //(
                    //IOobject
                    //(
                        //Foam::fvMesh::defaultRegion,
                        //runTime.timeName(),
                        //runTime,
                        //Foam::IOobject::MUST_READ
                    //)
                //)
            //);


        }

        ~SingleCellSamplerTest()
        {
            system("rm -r system");
            system("rm -r constant");
            system("rm -r 0");
        }
};

TEST_F(SingleCellSamplerTest, FullConstructor)
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
    SingleCellSampler sampler("SingleCellSampler", patch, 3.0);

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());

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
