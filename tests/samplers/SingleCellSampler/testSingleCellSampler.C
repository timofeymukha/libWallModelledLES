#include "fvCFD.H"
#include "NewtonRootFinder.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"

class SingleCellSamplerTest : public ::testing::Test
{
    protected:
        void SetUp() override
        {
            system("cp -r ../../testCases/channel_flow/system .");
            system("cp -r ../../testCases/channel_flow/constant .");
        }

        void TearDown() override
        {
            system("rm -r system");
            system("rm -r constant");
        }
};

TEST_F(SingleCellSamplerTest, FullConstructor)
{
    extern int mainArgc;
    extern char** mainArgv;
    int argc = mainArgc;
    char** argv = mainArgv;

//#include "setRootCase.H"
//Foam::argList::noBanner();
Foam::argList args(argc, argv);
Foam::Time runTime(Foam::Time::controlDictName, args);
Foam::polyMesh mesh
(
    Foam::IOobject
    (
	Foam::polyMesh::defaultRegion,
	runTime.timeName(),
	runTime,
	Foam::IOobject::MUST_READ
    )
);

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
