#include "codeRules.H"
#include "fvCFD.H"
#include "CaiSagautLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class CaiSagautLawOfTheWallTest : public ChannelFlow
{};

TEST_F(CaiSagautLawOfTheWallTest, ConstructFromConstants)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.395, 4);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 4);
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 4);
    CaiSagautLawOfTheWall law =
        CaiSagautLawOfTheWall("CaiSagautLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 4);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, ConstructFromTypeDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    CaiSagautLawOfTheWall law =
        CaiSagautLawOfTheWall("CaiSagautLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B(), 5.5);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 5.5);
}

TEST_F(CaiSagautLawOfTheWallTest, CopyConstructor)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.395, 4);
    CaiSagautLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, CopyAssignment)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.395, 4);
    CaiSagautLawOfTheWall law2  = CaiSagautLawOfTheWall(0.95, 2);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 4);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, Clone)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.395, 4);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<CaiSagautLawOfTheWall &>(law2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<CaiSagautLawOfTheWall &>(law2()).B(), 4);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 4);
}

TEST_F(CaiSagautLawOfTheWallTest, Value)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.41, 5.2);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, 0.01596792222649768);
}

TEST_F(CaiSagautLawOfTheWallTest, Derivative)
{
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.4, 5.5);

    scalar derivative = law.derivative(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, 1);
}

TEST_F(CaiSagautLawOfTheWallTest, ValueSampler)
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
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.41, 5.2);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, 0.01596792222649768);
}

TEST_F(CaiSagautLawOfTheWallTest, DerivativeSampler)
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
    CaiSagautLawOfTheWall law = CaiSagautLawOfTheWall(0.4, 5.5);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, 1);
}

