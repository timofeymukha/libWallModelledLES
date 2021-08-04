#include "codeRules.H"
#include "fvCFD.H"
#include "scalarListIOList.H"
#include "SampledVelocityField.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include <random>

class SampledVelocityTest : public ChannelFlow
{};

TEST_F(SampledVelocityTest, Constructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledVelocityField sampledField(patch);
}


TEST_F(SampledVelocityTest, Clone)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledVelocityField sampledField(patch);

    autoPtr<SampledField> clone(sampledField.clone());
}


TEST_F(SampledVelocityTest, NDims)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledVelocityField sampledField(patch);

    ASSERT_EQ(sampledField.nDims(), 3);
}


TEST_F(SampledVelocityTest, Name)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledVelocityField sampledField(patch);

    ASSERT_EQ(sampledField.name(), "U");
}


TEST_F(SampledVelocityTest, RegisterFieldsZero)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledVelocityField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("U"));

    const scalarListIOList & sampledFieldIOobject = sampledField.db().lookupObject<scalarListIOList>("U");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            ASSERT_EQ(sampledFieldIOobject[i][j], 0);
        }
    }
}


TEST_F(SampledVelocityTest, RegisterFieldsInitialize)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledVelocityField sampledField(patch, "cell");

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("U"));

    const scalarListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListIOList>("U");

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
                ASSERT_EQ
                (
                    sampledFieldIOobject[i][j],
                    U[patch.faceCells()[i]][j]
                );
            }
        }
    }
}


TEST_F(SampledVelocityTest, RegisterFieldsRead)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    system("cp -r 0/wallModelSamplingSingle 0/wallModelSampling");

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledVelocityField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("U"));

    const scalarListIOList & sampledFieldIOobject = sampledField.db().lookupObject<scalarListIOList>("U");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
                ASSERT_EQ(sampledFieldIOobject[i][j], i + 1);
        }
    }

}


TEST_F(SampledVelocityTest, Sample)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();


    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");
    // Init U to something varying and easily to test
    U.primitiveFieldRef() = mesh.C();

    SampledVelocityField sampledField(patch, "cellPointFace");

    labelList indexList(patch.faceCells());

    scalarField h(patch.size(), 0.19);


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
                ASSERT_FLOAT_EQ
                (
                    sampledValues[i][j],
                    U[indexList[i]][j]
                );
            }
        }
    }

}