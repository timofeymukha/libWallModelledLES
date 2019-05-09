#include "codeRules.H"
#include "fvCFD.H"
#include "SampledVelocityField.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"

class SampledVelocityTest : public ChannelFlow
{};

TEST_F(SampledVelocityTest, Constructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledVelocityField sampledField(patch);
}
