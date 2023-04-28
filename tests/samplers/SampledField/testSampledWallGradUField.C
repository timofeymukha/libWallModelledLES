#include "codeRules.H"
#include "fvCFD.H"
#include "scalarListIOList.H"
#include "SampledWallGradUField.H"
#undef Log
#include "gtest.h"
#include "fixtures.H"
#include <random>
#include "MultiCellSampler.H"

class SampledWallGradUTest : public ChannelFlow
{};

TEST_F(SampledWallGradUTest, Constructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledWallGradUField sampledField(patch);
}


TEST_F(SampledWallGradUTest, Clone)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledWallGradUField sampledField(patch);

    autoPtr<SampledField> clone(sampledField.clone());
}


TEST_F(SampledWallGradUTest, NDims)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledWallGradUField sampledField(patch);

    ASSERT_EQ(sampledField.nDims(), 3);
}


TEST_F(SampledWallGradUTest, Name)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledWallGradUField sampledField(patch);

    ASSERT_EQ(sampledField.name(), "wallGradU");
}


TEST_F(SampledWallGradUTest, RegisterFieldsZero)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);
    createSamplingHeightField(mesh);

    SampledWallGradUField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the global wallGradU field in the registry
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallGradU"));

    // Assert we registred the field in the wm registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("wallGradU"));

    const scalarListIOList & sampledFieldIOobject = sampledField.db().lookupObject<scalarListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            ASSERT_EQ(sampledFieldIOobject[i][j], 0);
        }
    }
}


TEST_F(SampledWallGradUTest, RegisterFieldsZeroMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    createSamplingHeightField(mesh);
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    h.boundaryFieldRef()[patch.index()] == 2;

    autoPtr<MultiCellSampler> sampler
    (
        new MultiCellSampler
        (
            "MultiCellSampler",
            patch,
            3.0,
            "cell",
            "Crawling",
            "CubeRootVol",
            true
        )
    );

    labelListList indexList = sampler->indexList();

    sampler->db().checkOut("wallGradU");
    sampler.clear();
    
    SampledWallGradUField sampledField(patch);
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("wallGradU"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i][0], j)
        {
            ASSERT_EQ(sampledFieldIOobject[i][0][j], 0);
        }
    }
}


TEST_F(SampledWallGradUTest, RegisterFieldsInitialize)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createVelocityField(mesh);
    createSamplingHeightField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U 
    U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledWallGradUField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the global wallGradU field in the registry
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallGradU"));

    // Assert we registred the field in wm the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("wallGradU"));

    const scalarListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            if (j == 1)
            {
                ASSERT_EQ(sampledFieldIOobject[i][j], 0);
            }
            else
            {
                ASSERT_NEAR
                (
                    sampledFieldIOobject[i][j],
                    U[patch.faceCells()[i]][j]/0.1,
                    1e-8
                );
            }
        }
    }
}


TEST_F(SampledWallGradUTest, RegisterInitializeMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    createSamplingHeightField(mesh);
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    h.boundaryFieldRef()[patch.index()] == 2;

    autoPtr<MultiCellSampler> sampler
    (
        new MultiCellSampler
        (
            "MultiCellSampler",
            patch,
            3.0,
            "cell",
            "Crawling",
            "CubeRootVol",
            true
        )
    );

    labelListList indexList = sampler->indexList();

    sampler->db().checkOut("wallGradU");
    sampler.clear();
    

    volVectorField & wallGradU = mesh.lookupObjectRef<volVectorField>("wallGradU");
    wallGradU.boundaryFieldRef()[patch.index()] == vector(1, 2, 3);

    SampledWallGradUField sampledField(patch);
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("wallGradU"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {

        ASSERT_EQ(sampledFieldIOobject[i].size(), 1);
        for(int j=0; j<sampledFieldIOobject[i][0].size(); j++)
        {
            
            if (j != 1)
            {
                ASSERT_FLOAT_EQ(sampledFieldIOobject[i][0][j], j+1);
            }
            else
            {
                ASSERT_FLOAT_EQ(sampledFieldIOobject[i][0][j], 0);
            }

        }
    }
}


TEST_F(SampledWallGradUTest, RegisterFieldsRead)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    system("cp -r 0/wallModelSamplingSingle 0/wallModelSampling");

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);
    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledWallGradUField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the global wallGradU field in the registry
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallGradU"));

    // Assert we registred the field in the wm registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("wallGradU"));

    const scalarListIOList & sampledFieldIOobject = sampledField.db().lookupObject<scalarListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
                ASSERT_EQ(sampledFieldIOobject[i][j], i + 1);
        }
    }

}


TEST_F(SampledWallGradUTest, RegisterReadFieldsMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    system("cp -r 0/wallModelSamplingMulti 0/wallModelSampling");

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    createSamplingHeightField(mesh);
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    h.boundaryFieldRef()[patch.index()] == 2;

    autoPtr<MultiCellSampler> sampler
    (
        new MultiCellSampler
        (
            "MultiCellSampler",
            patch,
            3.0,
            "cell",
            "Crawling",
            "CubeRootVol",
            true
        )
    );

    labelListList indexList = sampler->indexList();

    sampler->db().checkOut("wallGradU");
    sampler.clear();
    
    SampledWallGradUField sampledField(patch);
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("wallGradU"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("wallGradU");

    forAll(sampledFieldIOobject, i)
    {

        ASSERT_EQ(sampledFieldIOobject[i].size(), 1);
        for(int j=0; j<sampledFieldIOobject[i][0].size(); j++)
        {
            ASSERT_FLOAT_EQ(sampledFieldIOobject[i][0][j], j+1);
        }
    }
}
TEST_F(SampledWallGradUTest, Sample)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledWallGradUField sampledField(patch);

    // Init indexList with invalid ids since it is not used
    labelList indexList(patch.size(), -1);

    // Create wallGradU field
    const volScalarField & h = mesh.lookupObject<volScalarField>("hSampler");

    volVectorField & wallGradU =
        mesh.lookupObjectRef<volVectorField>("wallGradU");
    vectorField & boundaryValues =
        wallGradU.boundaryFieldRef()[patch.index()];
    boundaryValues = patch.Cf();

    scalarListList sampledValues(patch.size());

    sampledField.sample(sampledValues, indexList, h);


    forAll(sampledValues, i)
    {
        forAll(sampledValues[i], j)
        {
            if (j == 1)
            {
                ASSERT_EQ(sampledValues[i][j], 0);
            }
            else
            {
                ASSERT_NEAR
                (
                    sampledValues[i][j],
                    boundaryValues[i][j], 
                    1e-8
                );
            }
        }
    }

}


TEST_F(SampledWallGradUTest, SampleMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledWallGradUField sampledField(patch);

    // Init indexList with invalid ids since it is not used
    labelListList indexList(patch.size(), labelList(3, -1));

    volVectorField & wallGradU =
        mesh.lookupObjectRef<volVectorField>("wallGradU");
    vectorField & boundaryValues =
        wallGradU.boundaryFieldRef()[patch.index()];
    boundaryValues = patch.Cf();

    scalarListListList sampledValues(patch.size());

    sampledField.sample(sampledValues, indexList);

    forAll(sampledValues, i)
    {
        forAll(sampledValues[i][0], j)
        {
            if (j == 1)
            {
                ASSERT_FLOAT_EQ(sampledValues[i][0][j], 0);
            }
            else
            {
                ASSERT_NEAR
                (
                    sampledValues[i][0][j],
                    boundaryValues[i][j], 
                    1e-8
                );
            }
        }
    }

}


TEST_F(SampledWallGradUTest, Recompute)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();


    const fvPatch & patch = mesh.boundary()["topWall"];
    createWallModelSubregistry(mesh, patch);
    createVelocityField(mesh);
    createSamplingHeightField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    SampledWallGradUField sampledField(patch);


    // Init indexList with invalid ids since it is not used
    labelList indexList(patch.size(), -1);

    const volVectorField & wallGradU =
        mesh.lookupObject<volVectorField>("wallGradU");
    const vectorField & boundaryValues =
        wallGradU.boundaryField()[patch.index()];

    sampledField.recompute();


    forAll(boundaryValues, i)
    {
        forAll(boundaryValues[i], j)
        {
            ASSERT_NEAR
            (
                boundaryValues[i][j], 
                U[patch.faceCells()[i]][j]*patch.deltaCoeffs()[i],
                1e-8
            );
        }
    }

}
