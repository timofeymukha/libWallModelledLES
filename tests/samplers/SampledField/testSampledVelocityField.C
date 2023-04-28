#include "codeRules.H"
#include "fvCFD.H"
#include "scalarListIOList.H"
#include "SampledVelocityField.H"
#include "MultiCellSampler.H"
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


TEST_F(SampledVelocityTest, RegisterFieldsZeroMultiCell)
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
            "Crawling"
        )
    );

    labelListList indexList = sampler->indexList();

    // Remove U from the registry and delete the sampler
    sampler->db().checkOut("U");
    sampler.clear();
    

    h.boundaryFieldRef()[patch.index()] == 2;

    SampledVelocityField sampledField(patch);
    
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("U"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("U");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            forAll(sampledFieldIOobject[i][j], k)
            {
                ASSERT_EQ(sampledFieldIOobject[i][j][k], 0);
            }
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



TEST_F(SampledVelocityTest, RegisterFieldsInitializeMultiCell)
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

    createSamplingHeightField(mesh);
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );
    h.boundaryFieldRef()[patch.index()] == 3;

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

    // Remove U from the registry and delete the sampler
    sampler->db().checkOut("U");
    sampler.clear();
    
    h.boundaryFieldRef()[patch.index()] == 2;

    SampledVelocityField sampledField(patch);
    
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("U"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("U");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            forAll(sampledFieldIOobject[i][j], k)
            {
                if (k != 1)
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], U[indexList[i][j]][k]);
                }
                else
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], 0.0);
                }
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


TEST_F(SampledVelocityTest, RegisterFieldsReadMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    auto code = system("cp -r 0/wallModelSamplingMulti 0/wallModelSampling");
    ASSERT_EQ(code, 0);


    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    //U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];

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

    // Remove U from the registry and delete the sampler
    sampler->db().checkOut("U");
    sampler.clear();
    
    SampledVelocityField sampledField(patch);
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("U"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("U");

    // The mesh has 3 cells in x, 3 in z and 10 in y
    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            forAll(sampledFieldIOobject[i][j], k)
            {
                ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], k + 1);
            }
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

    SampledVelocityField sampledField(patch, "cell");

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


TEST_F(SampledVelocityTest, SampleMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];

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
            true,
            false
        )
    );

    labelListList indexList = sampler->indexList();

    // Remove U from the registry and delete the sampler
    sampler->db().checkOut("U");
    sampler.clear();
    
    SampledVelocityField sampledField(patch);

    createVelocityField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("U");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    
    sampledField.registerFields(indexList);

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("U");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            forAll(sampledFieldIOobject[i][j], k)
            {
                if (k != 1)
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], U[indexList[i][j]][k]);
                }
                else
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], 0.0);
                }
            }
        }
    }

}


TEST_F(SampledVelocityTest, CheckInterpolationWorks)
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
    
    forAll(U.primitiveFieldRef(), i)
    {
        for(int j=0; j<3; j++)
        {
            U.primitiveFieldRef()[i][j] = mesh.C()[i][1];
        }
    }

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
                    0.19
                );
            }
        }
    }

}


TEST_F(SampledVelocityTest, CheckInterpolationPointOutside)
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
    
    forAll(U.primitiveFieldRef(), i)
    {
        for(int j=0; j<3; j++)
        {
            U.primitiveFieldRef()[i][j] = mesh.C()[i][1];
        }
    }

    SampledVelocityField sampledField(patch, "cellPointFace");


    // we sample from near-wall cell which ends at 0.2
    // but set h to 0.5
    labelList indexList(patch.faceCells());
    scalarField h(patch.size(), 0.5);


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
                    0.2 // the top value of the cell we sample from
                );
            }
        }
    }

}