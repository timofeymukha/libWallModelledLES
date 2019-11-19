#include "codeRules.H"
#include "fvCFD.H"
#include "ReichardtLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "SingleCellSampler.H"

class ReichardtLawOfTheWallTest : public ChannelFlow
{};

TEST_F(ReichardtLawOfTheWallTest, ConstructFromConstants)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 7.8);

    dictionary dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);
}

TEST_F(ReichardtLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B1", 11);
    dict.add("B2", 3);
    dict.add("C", 7.8);


    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 7.8);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);
}  

TEST_F(ReichardtLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 7.8);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);
}  

TEST_F(ReichardtLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B1", 11);
    dict.add("B2", 3);
    dict.add("C", 7.8);
    ReichardtLawOfTheWall law =
        ReichardtLawOfTheWall("ReichardtLawOfTheWall", dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 7.8);

    dict = law.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);   
}

TEST_F(ReichardtLawOfTheWallTest, CopyConstructor)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    ReichardtLawOfTheWall law2(law);

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B1(), 11);
    ASSERT_DOUBLE_EQ(law2.B2(), 3);
    ASSERT_DOUBLE_EQ(law2.C(), 7.8);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);        
}

TEST_F(ReichardtLawOfTheWallTest, CopyAssignment)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    ReichardtLawOfTheWall law2 = ReichardtLawOfTheWall(1, 2, 4, 3);
    law2 = law;

    ASSERT_DOUBLE_EQ(law2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(law2.B1(), 11);
    ASSERT_DOUBLE_EQ(law2.B2(), 3);
    ASSERT_DOUBLE_EQ(law2.C(), 7.8);

    dictionary dict = law2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8); 
}

TEST_F(ReichardtLawOfTheWallTest, Clone)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<ReichardtLawOfTheWall &>(law2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<ReichardtLawOfTheWall &>(law2()).B1(), 11);
    ASSERT_DOUBLE_EQ(dynamic_cast<ReichardtLawOfTheWall &>(law2()).B2(), 3);
    ASSERT_DOUBLE_EQ(dynamic_cast<ReichardtLawOfTheWall &>(law2()).C(), 7.8);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);
}

TEST_F(ReichardtLawOfTheWallTest, Value)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar value = law.value(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -10.442820787986026);
}

TEST_F(ReichardtLawOfTheWallTest, Derivative)
{
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar derivative = law.derivative(0.5, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -375.6313131313131);
}

TEST_F(ReichardtLawOfTheWallTest, ValueSampler)
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
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -21.194402785132265);
}

TEST_F(ReichardtLawOfTheWallTest, DerivativeSampler)
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
    ReichardtLawOfTheWall law = ReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -62.97229219143577);
}

