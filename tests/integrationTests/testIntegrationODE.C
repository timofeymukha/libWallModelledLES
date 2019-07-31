#include "gtest.h"
#include <cstdlib>
#include <iostream>

class ChannelFlow : public ::testing::Test
{
    public:
        ChannelFlow()
        {
            std::system("cp -r ../testCases/channel_flow/system .");
            std::system("cp -r ../testCases/channel_flow/constant .");
            std::system("cp -r ../testCases/channel_flow/0 .");
        }

        virtual ~ChannelFlow()
        {
            std::system("rm -r system");
            std::system("rm -r constant");
            std::system("rm -r 0");
            std::system("rm -r 0.01");
            std::system("rm -r processor*");
        }
};

class IntegrationTest : public ChannelFlow
{};

// * * * * * * * * * * * * * * * Run LOTW in serial * * * * * * * * * * * * * //

TEST_F(IntegrationTest, RunEquilibriumODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutEquilibriumODEVanDriest");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunEquilibriumODEDurat)
{
    int success = std::system("changeDictionary -dict system/setNutEquilibriumODEDuprat");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunPGradODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutPGradODEVanDriest");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunPGradODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutPGradODEDuprat");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * decomposePar with ODE * * * * * * * * * * * * * //

TEST_F(IntegrationTest, DecomposeEquilibriumODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutEquilibriumODEVanDriest");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeEquilibriumODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutEquilibriumODEDuprat");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposePGradODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutPGradODEDuprat");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposePGradODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutPGradODEVanDriest");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}


// * * * * * * * * * * * * * * * Run in parallel with ODE * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ParallelRunEquilibuimODEVanDriest)
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
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunEquilibuimODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutEquilibriumODEDuprat -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunPGradODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutPGradODEDuprat -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunPGradODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutPGradODEVanDriest -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * Reconstruct with ODE * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ReconstructEquilibriumODEVanDriest)
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

TEST_F(IntegrationTest, ReconstructEquilibriumODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutEquilibriumODEDuprat -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructPGradODEVanDriest)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutPGradODEVanDriest -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructPGradODEDuprat)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutPGradODEDuprat -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}