#include "gtest.h"
#include <cstdlib>
#include <iostream>

class ChannelFlow : public ::testing::Test
{
    public:
        ChannelFlow()
        {
            auto code = std::system("cp -r ../testCases/channel_flow/system .");
            code = std::system("cp -r ../testCases/channel_flow/constant .");
            code = std::system("cp -r ../testCases/channel_flow/0 .");
        }

        virtual ~ChannelFlow()
        {
            auto code = std::system("rm -r system");
            code = std::system("rm -r constant");
            code = std::system("rm -r 0");
            code = std::system("rm -r 0.01");
            code = std::system("rm -r processor*");
        }
};

class IntegrationTest : public ChannelFlow
{};

// * * * * * * * * * * * * * * * Run in serial * * * * * * * * * * * * * //

TEST_F(IntegrationTest, RunKnownWallShearStress)
{
    int success = std::system("changeDictionary -dict system/setNutKnownWallShearStress");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * decomposePar * * * * * * * * * * * * * //

TEST_F(IntegrationTest, DecomposeKnownWallShearStress)
{
    int success = std::system("changeDictionary -dict system/setNutKnownWallShearStress");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * Run in parallel * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ParallelRunKnowWallShearStress)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutKnownWallShearStress -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * Reconstruct * * * * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ReconstructKnownWallShearStress)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutEquilibriumODEVanDriest -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}