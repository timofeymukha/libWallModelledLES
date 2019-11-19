#include "codeRules.H"
#include "fvCFD.H"
#include "SpaldingLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class SpaldingLawOfTheWallTest : public ChannelFlow
{};

TEST_F(SpaldingLawOfTheWallTest, ConstructFromConstants)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.395, 4);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(SpaldingLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 4);
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);    
}

TEST_F(SpaldingLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 4);
    SpaldingLawOfTheWall law =
        SpaldingLawOfTheWall("SpaldingLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);     
}

TEST_F(SpaldingLawOfTheWallTest, ConstructFromTypeDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    SpaldingLawOfTheWall law =
        SpaldingLawOfTheWall("SpaldingLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B(), 5.5);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 5.5);     
}

TEST_F(SpaldingLawOfTheWallTest, CopyConstructor)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.395, 4);
    SpaldingLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);    
}

TEST_F(SpaldingLawOfTheWallTest, CopyAssignment)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.395, 4);
    SpaldingLawOfTheWall law2  = SpaldingLawOfTheWall(0.95, 2);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(SpaldingLawOfTheWallTest, Clone)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.395, 4);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<SpaldingLawOfTheWall &>(law2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<SpaldingLawOfTheWall &>(law2()).B(), 4);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);   
}

TEST_F(SpaldingLawOfTheWallTest, Value)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.4, 5.5);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -975.4136107911548);
}

TEST_F(SpaldingLawOfTheWallTest, Derivative)
{
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.4, 5.5);

    scalar derivative = law.derivative(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -27111.848542674237);
}

TEST_F(SpaldingLawOfTheWallTest, ValueSampler)
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
        "SingleCellSampler",
        patch,
        3.0
    );
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.4, 5.5);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -500);
}

TEST_F(SpaldingLawOfTheWallTest, DerivativeSampler)
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
        "SingleCellSampler",
        patch,
        3.0
    );
    SpaldingLawOfTheWall law = SpaldingLawOfTheWall(0.4, 5.5);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -12500);
}

