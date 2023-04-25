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

// * * * * * * * * * * * * * * * Run LOTW in serial * * * * * * * * * * * * * //

TEST_F(IntegrationTest, RunLOTWSpalding)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWSpalding");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunLOTWReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWReichardt");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunLOTWIntegratedReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWIntegratedReichardt");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunLOTWWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWWernerWengle");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunLOTWRoughLogLaw)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWRoughLogLaw");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, RunLOTWIntegratedWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWWernerWengle");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("pimpleFoam");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * decomposePar with LOTW * * * * * * * * * * * * * //

TEST_F(IntegrationTest, DecomposeLOTWSpalding)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWSpalding");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeLOTWReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWReichardt");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeLOTWIntegratedReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWIntegratedReichardt");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeLOTWWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWWernerWengle");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeLOTWIntegratedWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWWernerWengle");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, DecomposeLOTWIntegratedRoughLogLaw)
{
    int success = std::system("changeDictionary -dict system/setNutLOTWRoughLogLaw");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * Run in parallel with LOTW * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ParallelRunLOTWSpalding)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWSpalding -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunLOTWReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWReichardt -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunLOTWWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWWernerWengle -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunLOTWIntegratedWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWIntegratedWernerWengle -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunLOTWIntegratedReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWIntegratedReichardt -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ParallelRunLOTWRoughLogLaw)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWRoughLogLaw -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 pimpleFoam -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

// * * * * * * * * * * * * * * * Reconstruct with LOTW * * * * * * * * * * * * * //

TEST_F(IntegrationTest, ReconstructLOTWSpalding)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWSpalding -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructLOTWReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWReichardt -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructLOTWIntegratedReichardt)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWIntegratedReichardt -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructLOTWWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWWernerWengle -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructLOTWIntegratedWernerWengle)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWIntegratedWernerWengle -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}

TEST_F(IntegrationTest, ReconstructLOTWRoughLogLaw)
{
    int success = std::system("changeDictionary -dict system/setNutFixedValue");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("decomposePar -force");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("mpirun -np 2 changeDictionary -dict system/setNutLOTWRoughLogLaw -parallel");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
    success = std::system("reconstructPar -withZero");
    ASSERT_EQ(WIFEXITED(success), true);
    ASSERT_EQ(WEXITSTATUS(success), 0);
}