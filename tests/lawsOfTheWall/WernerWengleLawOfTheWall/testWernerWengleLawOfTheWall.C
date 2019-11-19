#include "codeRules.H"
#include "fvCFD.H"
#include "WernerWengleLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class WernerWengleLawOfTheWallTest : public ChannelFlow
{};

TEST_F(WernerWengleLawOfTheWallTest, ConstructFromConstants)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(0.395, 4);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(WernerWengleLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("A", 0.395);
    dict.add("B", 4);
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);    
}

TEST_F(WernerWengleLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.A(), 8.3);
    ASSERT_FLOAT_EQ(law.B(), 0.14285714285714285);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 8.3);
    ASSERT_FLOAT_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 0.142857);    
}

TEST_F(WernerWengleLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("A", 0.395);
    dict.add("B", 4);
    WernerWengleLawOfTheWall law =
        WernerWengleLawOfTheWall("WernerWengleLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.A(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);     
}

TEST_F(WernerWengleLawOfTheWallTest, CopyConstructor)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(0.395, 4);
    WernerWengleLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.A(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);    
}

TEST_F(WernerWengleLawOfTheWallTest, CopyAssignment)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(0.395, 4);
    WernerWengleLawOfTheWall law2  = WernerWengleLawOfTheWall(0.95, 2);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.A(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(WernerWengleLawOfTheWallTest, Clone)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(0.395, 4);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<WernerWengleLawOfTheWall &>(law2()).A(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<WernerWengleLawOfTheWall &>(law2()).B(), 4);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("A", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);   
}

TEST_F(WernerWengleLawOfTheWallTest, Value)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(8.3, 1./7);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -9.766375100821723);
}

TEST_F(WernerWengleLawOfTheWallTest, Derivative)
{
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(8.3, 1./7);

    scalar derivative = law.derivative(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -392.02276821722046);
}

TEST_F(WernerWengleLawOfTheWallTest, ValueSampler)
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
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(8.3, 1./7);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -20.167182846190865);
}

TEST_F(WernerWengleLawOfTheWallTest, DerivativeSampler)
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
    WernerWengleLawOfTheWall law = WernerWengleLawOfTheWall(8.3, 1./7);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -72.02565302211023);
}

