#include "codeRules.H"
#include "fvCFD.H"
#include "IntegratedWernerWengleLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class IntegratedWernerWengleLawOfTheWallTest : public ChannelFlow
{};

TEST_F(IntegratedWernerWengleLawOfTheWallTest, ConstructFromConstants)
{
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(0.395, 4);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("A", 0.395);
    dict.add("B", 4);
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.A(), 8.3);
    ASSERT_FLOAT_EQ(law.B(), 0.14285714285714285);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 8.3);
    ASSERT_FLOAT_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 0.142857);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("A", 0.395);
    dict.add("B", 4);
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall("IntegratedWernerWengleLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, CopyConstructor)
{
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(0.395, 4);
    IntegratedWernerWengleLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.A(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, CopyAssignment)
{
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(0.395, 4);
    IntegratedWernerWengleLawOfTheWall law2 =
        IntegratedWernerWengleLawOfTheWall(0.95, 2);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.A(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, Clone)
{
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(0.395, 4);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedWernerWengleLawOfTheWall &>(law2()).A(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedWernerWengleLawOfTheWall &>(law2()).B(), 4);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, Value)
{
    IntegratedWernerWengleLawOfTheWall law = 
        IntegratedWernerWengleLawOfTheWall(8.3, 1./7);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, 0.01277344612632577);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, Derivative)
{
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(8.3, 1./7);

    scalar derivative = law.derivative();
    ASSERT_DOUBLE_EQ(derivative, 1);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, ValueSampler)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        patch,
        3.0
    );
    IntegratedWernerWengleLawOfTheWall law =
        IntegratedWernerWengleLawOfTheWall(8.3, 1./7);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_FLOAT_EQ(value, 0.039920206660171285);
}

TEST_F(IntegratedWernerWengleLawOfTheWallTest, DerivativeSampler)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        patch,
        3.0
    );
    IntegratedWernerWengleLawOfTheWall law = 
        IntegratedWernerWengleLawOfTheWall(8.3, 1./7);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, 1);
}

