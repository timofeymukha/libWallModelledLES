#include "codeRules.H"
#include "fvCFD.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"


// construct h given mesh
void createSamplingHeightField(const Foam::fvMesh & mesh);

// construct U given mesh
void createVelocityField(const Foam::fvMesh & mesh);

// Create a fvMesh
Foam::autoPtr<Foam::fvMesh> createMesh(const Foam::Time & runTime);

// Create subregistry for sampled fields for a given patch
void createWallModelSubregistry(const fvMesh & mesh, const fvPatch & patch);

class ChannelFlow : public ::testing::Test
{


    public:
        ChannelFlow()
        {
            system("cp -r testCases/channel_flow/system .");
            system("cp -r testCases/channel_flow/constant .");
            system("cp -r testCases/channel_flow/0 .");

        }

        ~ChannelFlow()
        {
            system("rm -r system");
            system("rm -r constant");
            system("rm -r 0");
        }
};

