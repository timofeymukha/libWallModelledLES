#include "codeRules.H"
#include "fvCFD.H"
#include "DupratEddyViscosity.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include "scalarListIOList.H"

class DupratEddyViscosityTest : public ChannelFlow
{};

TEST_F(DupratEddyViscosityTest, ConstructFromConstants)
{
    DupratEddyViscosity eddy = DupratEddyViscosity(0.395, 4, 0.78);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);
    ASSERT_DOUBLE_EQ(eddy.beta(), 0.78);

    dictionary dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 0.78);    
}

TEST_F(DupratEddyViscosityTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("APlus", 4);
    dict.add("beta", 0.78);    
    DupratEddyViscosity eddy = DupratEddyViscosity(dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);
    ASSERT_DOUBLE_EQ(eddy.beta(), 0.78);

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 0.78);    
}

TEST_F(DupratEddyViscosityTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    DupratEddyViscosity eddy = DupratEddyViscosity(dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 18);
    ASSERT_DOUBLE_EQ(eddy.beta(), 0.78);    

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 18);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 0.78);
}

TEST_F(DupratEddyViscosityTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("APlus", 4);
    dict.add("beta", 5);
    DupratEddyViscosity eddy =
        DupratEddyViscosity("DupratEddyViscosity", dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);
    ASSERT_DOUBLE_EQ(eddy.beta(), 5);

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 5);
}

TEST_F(DupratEddyViscosityTest, CopyConstructor)
{
    DupratEddyViscosity eddy = DupratEddyViscosity(0.395, 4, 0.5);
    DupratEddyViscosity eddy2(eddy);

    ASSERT_DOUBLE_EQ(eddy2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy2.APlus(), 4);
    ASSERT_DOUBLE_EQ(eddy2.beta(), 0.5);

    dictionary dict = eddy2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 0.5);
}

TEST_F(DupratEddyViscosityTest, CopyAssignment)
{
    DupratEddyViscosity eddy = DupratEddyViscosity(0.395, 4, 0.5);
    DupratEddyViscosity eddy2  = DupratEddyViscosity(0.95, 20, 0.8);
    eddy2 = eddy;

    ASSERT_DOUBLE_EQ(eddy2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy2.APlus(), 4);
    ASSERT_DOUBLE_EQ(eddy2.beta(), 0.5);

    dictionary dict = eddy2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 0.5);
}

TEST_F(DupratEddyViscosityTest, Clone)
{
    DupratEddyViscosity eddy = DupratEddyViscosity(0.395, 4, 1);
    autoPtr<EddyViscosity> eddy2  = eddy.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<DupratEddyViscosity &>(eddy2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<DupratEddyViscosity &>(eddy2()).APlus(), 4);
    ASSERT_DOUBLE_EQ(dynamic_cast<DupratEddyViscosity &>(eddy2()).beta(), 1);

    dictionary dict = eddy2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("beta", 0.0), 1);
}

TEST_F(DupratEddyViscosityTest, AddFieldsToSampler)
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
        3.0,
        "cell",
        "Tree", 
        "CubeRootVol",
        false
    );

    DupratEddyViscosity eddy = DupratEddyViscosity(0.4, 18, 0.78);
    eddy.addFieldsToSampler(sampler);

    ASSERT_EQ(sampler.nSampledFields(), 3);
    ASSERT_TRUE(sampler.db().foundObject<scalarListIOList>("pGrad"));

}

TEST_F(DupratEddyViscosityTest, Value)
{
    DupratEddyViscosity eddy = DupratEddyViscosity(dictionary());
    scalarList y(2, 0.01);
    y[1] = 0.1;

    scalarList values = eddy.value(y, 0.1, 0.04, 8e-6);
    scalarList refValues(2);
    refValues[0] = 0.00021063538749341106;
    refValues[1] = 0.007392720890351022;

    ASSERT_DOUBLE_EQ(values[0], refValues[0]);
    ASSERT_DOUBLE_EQ(values[1], refValues[1]);
}

TEST_F(DupratEddyViscosityTest, ValueSampler)
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
    DupratEddyViscosity eddy = DupratEddyViscosity(0.4, 18, 0.78);
    eddy.addFieldsToSampler(sampler);

    scalarList y(2, 0.01);
    y[1] = 0.1;
    scalarList values = eddy.value(sampler, 5, y, 0.04, 8e-6);

    scalarList refValues(2);
    refValues[0] = 0.00013779990983794432;
    refValues[1] = 0.0015999999999880782;

    ASSERT_DOUBLE_EQ(values[0], refValues[0]);
    ASSERT_DOUBLE_EQ(values[1], refValues[1]);
}

