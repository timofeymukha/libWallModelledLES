#include "codeRules.H"
#include "fvCFD.H"
#include "IntegratedReichardtLawOfTheWall.H"
#include "SingleCellSampler.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"

class IntegratedReichardtLawOfTheWallTest : public ChannelFlow
{};

TEST_F(IntegratedReichardtLawOfTheWallTest, ConstructFromConstants)
{
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);

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

TEST_F(IntegratedReichardtLawOfTheWallTest, ConstructFromDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B1", 11);
    dict.add("B2", 3);
    dict.add("C", 7.8);


    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(dict);

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

TEST_F(IntegratedReichardtLawOfTheWallTest, ConstructFromDictionaryDefaultValues)
{
    dictionary dict = dictionary();

    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(dict);

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

TEST_F(IntegratedReichardtLawOfTheWallTest, ConstructFromTypeDictionary)
{
    dictionary dict = dictionary();
    dict.add("kappa", 0.395);
    dict.add("B1", 11);
    dict.add("B2", 3);
    dict.add("C", 7.8);
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall("IntegratedReichardtLawOfTheWall", dict);

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

TEST_F(IntegratedReichardtLawOfTheWallTest, CopyConstructor)
{
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    IntegratedReichardtLawOfTheWall law2(law);

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

TEST_F(IntegratedReichardtLawOfTheWallTest, CopyAssignment)
{
    IntegratedReichardtLawOfTheWall law = 
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    IntegratedReichardtLawOfTheWall law2 =
        IntegratedReichardtLawOfTheWall(1, 2, 4, 3);
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

TEST_F(IntegratedReichardtLawOfTheWallTest, Clone)
{
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);
    autoPtr<LawOfTheWall> law2  = law.clone();

    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedReichardtLawOfTheWall &>(law2()).kappa(), 0.395);
    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedReichardtLawOfTheWall &>(law2()).B1(), 11);
    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedReichardtLawOfTheWall &>(law2()).B2(), 3);
    ASSERT_DOUBLE_EQ(dynamic_cast<IntegratedReichardtLawOfTheWall &>(law2()).C(), 7.8);

    dictionary dict = law2().constDict();
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("kappa", 0.0), 0.395);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B1", 0.0), 11);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("B2", 0.0), 3);
    ASSERT_DOUBLE_EQ(dict.lookupOrDefault<scalar>("C", 0.0), 7.8);
}

TEST_F(IntegratedReichardtLawOfTheWallTest, Value)
{
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar value = law.value(0.5, 0.1, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(value, -0.03867378388240429);
}

TEST_F(IntegratedReichardtLawOfTheWallTest, Derivative)
{
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar derivative = law.derivative(0.1, 0.2, 0.04, 8e-6);
    ASSERT_DOUBLE_EQ(derivative, -2.4691238790839787);
}

TEST_F(IntegratedReichardtLawOfTheWallTest, ValueSampler)
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
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar value = law.value(sampler, 5, 0.04, 8e-6);
    ASSERT_FLOAT_EQ(value, -0.20040621277606402);
}


TEST_F(IntegratedReichardtLawOfTheWallTest, ValueMulticellSampler)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    for (int i=0; i< U.size(); i++)
    {
        U.primitiveFieldRef()[i] = vector(5, 0, 0);
    }

    h.boundaryFieldRef()[patch.index()] == 2;
    MultiCellSampler sampler
    (
        patch,
        3.0,
        "cell",
        "Crawling",
        "WallNormalDistance",
        true,
        false
    );
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.4, 11, 3, 7.8);

    // label index = 5;
    // const scalarListList & sampledU =
    //     sampler.db().lookupObject<scalarListListIOList>("U")[index];
    // scalarList samplerh = sampler.h()[index];
    // scalarList samplerl = sampler.lengthList()[index];
    
    scalar value = law.valueMulticell(sampler, 5, 0.04, 8e-6);
    ASSERT_FLOAT_EQ(value, 1.6481687);
}

TEST_F(IntegratedReichardtLawOfTheWallTest, DerivativeSampler)
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
    IntegratedReichardtLawOfTheWall law =
        IntegratedReichardtLawOfTheWall(0.395, 11, 3, 7.8);

    scalar derivative = law.derivative(sampler, 5, 0.04, 8e-6);
    ASSERT_FLOAT_EQ(derivative, -5.515923941988172);
}

