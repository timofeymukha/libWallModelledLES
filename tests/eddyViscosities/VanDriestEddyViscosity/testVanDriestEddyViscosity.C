#include "codeRules.H"
#include "fvCFD.H"
#include "VanDriestEddyViscosity.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"

class VanDriestEddyViscosityTest : public ChannelFlow
{};

TEST_F(VanDriestEddyViscosityTest, ConstructFromConstants)
{
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(0.395, 4);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);

    dictionary dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
}

TEST_F(VanDriestEddyViscosityTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("APlus", 4);
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);    
}

TEST_F(VanDriestEddyViscosityTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 18);

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.4);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 18);    
}

TEST_F(VanDriestEddyViscosityTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("APlus", 4);
    VanDriestEddyViscosity eddy =
        VanDriestEddyViscosity("VanDriestEddyViscosity", dict);

    ASSERT_DOUBLE_EQ(eddy.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy.APlus(), 4);

    dict = eddy.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);     
}

TEST_F(VanDriestEddyViscosityTest, CopyConstructor)
{
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(0.395, 4);
    VanDriestEddyViscosity eddy2(eddy);

    ASSERT_DOUBLE_EQ(eddy2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy2.APlus(), 4);

    dictionary dict = eddy2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);    
}

TEST_F(VanDriestEddyViscosityTest, CopyAssignment)
{
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(0.395, 4);
    VanDriestEddyViscosity eddy2  = VanDriestEddyViscosity(0.95, 2);
    eddy2 = eddy;

    ASSERT_DOUBLE_EQ(eddy2.kappa(), 0.395);
    ASSERT_DOUBLE_EQ(eddy2.APlus(), 4);

    dictionary dict = eddy2.constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);
}

TEST_F(VanDriestEddyViscosityTest, Clone)
{
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(0.395, 4);
    autoPtr<EddyViscosity> eddy2  = eddy.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<VanDriestEddyViscosity &>(eddy2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<VanDriestEddyViscosity &>(eddy2()).APlus(), 4);

    dictionary dict = eddy2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("APlus", 0.0), 4);   
}

TEST_F(VanDriestEddyViscosityTest, Value)
{
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(dictionary());
    scalarList y(2, 0.01);
    y[1] = 0.1;

    scalarList values = eddy.value(y, 0.04, 8e-6);
    scalarList refValues(2);
    refValues[0] = 0.00014072205953523847;
    refValues[1] = 0.0015999999999972369;

    ASSERT_DOUBLE_EQ(values[0], refValues[0]);
    ASSERT_DOUBLE_EQ(values[1], refValues[1]);
}

TEST_F(VanDriestEddyViscosityTest, ValueSampler)
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
    VanDriestEddyViscosity eddy = VanDriestEddyViscosity(0.4, 18);
    
    scalarList y(2, 0.01);
    y[1] = 0.1;
    scalarList values = eddy.value(sampler, 5, y, 0.04, 8e-6);

    scalarList refValues(2);
    refValues[0] = 0.00014072205953523847;
    refValues[1] = 0.0015999999999972369;

    ASSERT_DOUBLE_EQ(values[0], refValues[0]);
    ASSERT_DOUBLE_EQ(values[1], refValues[1]);
}

