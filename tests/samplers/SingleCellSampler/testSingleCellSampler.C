#include "fvCFD.H"
#include "SingleCellSampler.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"

class SingleCellSamplerTest : public ::testing::Test
{
    public:
        Foam::autoPtr<Foam::fvMesh> mesh_;

        SingleCellSamplerTest()
        {
            //extern int mainArgc;
            //extern char** mainArgv;
            //int argc = mainArgc;
            //char** argv = mainArgv;

            //system("cp -r ../../testCases/channel_flow/system .");
            //system("cp -r ../../testCases/channel_flow/constant .");
            //system("cp -r ../../testCases/channel_flow/0 .");

            Foam::argList::noBanner();
            //Foam::argList args(argc, argv);
            //Foam::Time runTime(Foam::Time::controlDictName, args);
            //mesh_.reset
            //(
                //new Foam::fvMesh
                //(
                    //Foam::IOobject
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
            Info << "helj" <<nl;
        }

        void TearDown() override
        {
            Info << "tear" << nl;
        }
};

TEST_F(SingleCellSamplerTest, FullConstructor)
{
            extern int mainArgc;
            extern char** mainArgv;
            int argc = mainArgc;
            char** argv = mainArgv;
#include "setRootCase.H"
            Info << "hi" << nl;
    //const fvPatch & patch = mesh_->boundary()["bottomWall"];
    //SingleCellSampler sampler("SingleCellSampler", patch, 3.0);

    //ASSERT_EQ(sampler.averagingTime(), 3.0);
    //ASSERT_EQ(&sampler.patch(), &patch);
    //ASSERT_EQ(sampler.indexList().size(), patch.size());
    //ASSERT_EQ(sampler.lengthList().size(), patch.size());
    //ASSERT_EQ(sampler.h().size(), patch.size());
}

int mainArgc;
char** mainArgv;

int main(int argc, char **argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    ::testing::InitGoogleTest(&argc, argv);

    mainArgc = argc;
    mainArgv = argv;

    return RUN_ALL_TESTS();
}
