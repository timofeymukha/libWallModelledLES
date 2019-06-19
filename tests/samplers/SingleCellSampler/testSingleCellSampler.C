#include "fvCFD.H"
#include "SingleCellSampler.H"
#include <functional>
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"
#include "fixtures.H"


class SingleCellSamplerTest : public ChannelFlow
{};

TEST_F(SingleCellSamplerTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler("SingleCellSampler", patch, 3.0);

    ASSERT_EQ(&sampler.patch(), &patch);
    ASSERT_EQ(sampler.averagingTime(), 3.0);
    ASSERT_EQ(&sampler.mesh(), &mesh);
    ASSERT_EQ(sampler.indexList().size(), patch.size());
    ASSERT_EQ(sampler.lengthList().size(), patch.size());
    ASSERT_EQ(sampler.h().size(), patch.size());
}
