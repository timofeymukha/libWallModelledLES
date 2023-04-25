#include "codeRules.H"
#include "fvCFD.H"
#include "RoughLogLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class RoughLogLawOfTheWallTest : public ChannelFlow
{};

TEST_F(RoughLogLawOfTheWallTest, ConstructFromConstants)
{
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 11);
    ASSERT_DOUBLE_EQ(law.ks(), 3);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}

TEST_F(RoughLogLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 11);
    dict.add("ks", 3);


    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 11);
    ASSERT_DOUBLE_EQ(law.ks(), 3);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}  

TEST_F(RoughLogLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    dict.add("B", 11);
    dict.add("ks", 3);
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B(), 11);
    ASSERT_DOUBLE_EQ(law.ks(), 3);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}  

TEST_F(RoughLogLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B", 11);
    dict.add("ks", 3);
    RoughLogLawOfTheWall law =
        RoughLogLawOfTheWall("RoughLogLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B(), 11);
    ASSERT_DOUBLE_EQ(law.ks(), 3);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}

TEST_F(RoughLogLawOfTheWallTest, CopyConstructor)
{
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);
    RoughLogLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 11);
    ASSERT_DOUBLE_EQ(law2.ks(), 3);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}

TEST_F(RoughLogLawOfTheWallTest, CopyAssignment)
{
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);
    RoughLogLawOfTheWall law2 = RoughLogLawOfTheWall(1, 2, 4);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B(), 11);
    ASSERT_DOUBLE_EQ(law2.ks(), 3);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}

TEST_F(RoughLogLawOfTheWallTest, Clone)
{
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<RoughLogLawOfTheWall &>(law2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<RoughLogLawOfTheWall &>(law2()).B(), 11);
    ASSERT_DOUBLE_EQ(dynamic_cast<RoughLogLawOfTheWall &>(law2()).ks(), 3);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("ks", 0.0), 3);
}

TEST_F(RoughLogLawOfTheWallTest, Value)
{
    scalar u = 1;
    scalar kappa = 0.4;
    scalar B = 5;
    scalar ks = 0.1;
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(kappa, B, ks);

    scalar value = law.value(u, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, 0.04 - u/(1/kappa*Foam::log(0.2/ks) + B));
}

TEST_F(RoughLogLawOfTheWallTest, Derivative)
{
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);

    scalar derivative = law.derivative();
    ASSERT_DOUBLE_EQ(derivative, 1);
}

TEST_F(RoughLogLawOfTheWallTest, ValueSampler)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    h.boundaryFieldRef()[patch.index()] == 0;
    auto & U = const_cast<volVectorField &>
    (
        mesh.thisDb().lookupObject<volVectorField>("U")
    );

    U.primitiveFieldRef() = vector(1, 0, 0);

    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0
    );
    scalar u = 1;
    scalar kappa = 0.4;
    scalar B = 5;
    scalar ks = 0.1;
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(kappa, B, ks);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, 0.04 - u/(1/kappa*Foam::log(0.1/ks) + B));
}

TEST_F(RoughLogLawOfTheWallTest, DerivativeSampler)
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
    RoughLogLawOfTheWall law = RoughLogLawOfTheWall(0.395, 11, 3);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, 1);
}

