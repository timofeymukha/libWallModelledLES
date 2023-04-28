#include "codeRules.H"
#include "fvCFD.H"
#include "scalarListIOList.H"
#include "SampledPGradField.H"
#include "MultiCellSampler.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"
#include <random>

class SampledPGradTest : public ChannelFlow
{};

TEST_F(SampledPGradTest, Constructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledPGradField sampledField(patch);
}


TEST_F(SampledPGradTest, Clone)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledPGradField sampledField(patch);

    autoPtr<SampledField> clone(sampledField.clone());
}


TEST_F(SampledPGradTest, NDims)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledPGradField sampledField(patch);

    ASSERT_EQ(sampledField.nDims(), 3);
}


TEST_F(SampledPGradTest, Name)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SampledPGradField sampledField(patch);

    ASSERT_EQ(sampledField.name(), "pGrad");
}


TEST_F(SampledPGradTest, RegisterFieldsZero)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);
    createSamplingHeightField(mesh);

    SampledPGradField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("pGrad"));

    const scalarListIOList & sampledFieldIOobject = 
        sampledField.db().lookupObject<scalarListIOList>("pGrad");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            ASSERT_EQ(sampledFieldIOobject[i][j], 0);
        }
    }
}


TEST_F(SampledPGradTest, RegisterFieldsZeroMultiCell)
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

    // Remove U from the registry and delete the sampler
    sampler.clear();
    

    SampledPGradField sampledField(patch);
    
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("pGrad"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("pGrad");

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


TEST_F(SampledPGradTest, RegisterFieldsInitialize)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);
    createPGradField(mesh);
    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");

    pGrad.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledPGradField sampledField(patch, "cell");

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("pGrad"));

    const scalarListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListIOList>("pGrad");

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
                    pGrad[patch.faceCells()[i]][j]
                );
            }
        }
    }
}



TEST_F(SampledPGradTest, RegisterFieldsInitializeMultiCell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);
    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("hSampler")
    );

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    h.boundaryFieldRef()[patch.index()] == 3;

    createPGradField(mesh);
    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");

    pGrad.primitiveFieldRef() = mesh.C();

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
    sampler.clear();
    
    h.boundaryFieldRef()[patch.index()] == 3;

    SampledPGradField sampledField(patch);
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("pGrad"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("pGrad");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
            forAll(sampledFieldIOobject[i][j], k)
            {
                if (k != 1)
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k],
                                    pGrad[indexList[i][j]][k]);
                }
                else
                {
                    ASSERT_FLOAT_EQ(sampledFieldIOobject[i][j][k], 0.0);
                }
            }
        }
    }
}


TEST_F(SampledPGradTest, RegisterFieldsRead)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    system("cp -r 0/wallModelSamplingSingle 0/wallModelSampling");

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);
    createPGradField(mesh);
    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");

    pGrad.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);

    SampledPGradField sampledField(patch);

    sampledField.registerFields(patch.faceCells());

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListIOList>("pGrad"));

    const scalarListIOList & sampledFieldIOobject = sampledField.db().lookupObject<scalarListIOList>("pGrad");

    forAll(sampledFieldIOobject, i)
    {
        forAll(sampledFieldIOobject[i], j)
        {
                ASSERT_EQ(sampledFieldIOobject[i][j], i + 1);
        }
    }

}


TEST_F(SampledPGradTest, RegisterFieldsReadMulticell)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    // Make previously sampled data readable
    auto code = system("cp -r 0/wallModelSamplingMulti 0/wallModelSampling");
    ASSERT_EQ(code, 0);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    createSamplingHeightField(mesh);
    createPGradField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("pGrad");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    const fvPatch & patch = mesh.boundary()["bottomWall"];

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
    sampler->db().checkOut("pGrad");
    sampler.clear();
    
    SampledPGradField sampledField(patch);
    
    sampledField.registerFields(indexList);

    // Assert we registred the field in the registry
    ASSERT_TRUE(sampledField.db().foundObject<scalarListListIOList>("pGrad"));

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("pGrad");

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


TEST_F(SampledPGradTest, Sample)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);
    createSamplingHeightField(mesh);

    createPGradField(mesh);
    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");
    pGrad.primitiveFieldRef() = mesh.C();

    SampledPGradField sampledField(patch, "cell");

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
                    pGrad[indexList[i]][j]
                );
            }
        }
    }

}


TEST_F(SampledPGradTest, SampleMulticell)
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
    sampler->db().checkOut("pGrad");
    sampler.clear();
    
    SampledPGradField sampledField(patch);

    createPGradField(mesh);
    volVectorField & U = mesh.lookupObjectRef<volVectorField>("pGrad");

    // Init U to something varying
    U.primitiveFieldRef() = mesh.C();

    
    sampledField.registerFields(indexList);

    const scalarListListIOList & sampledFieldIOobject =
        sampledField.db().lookupObject<scalarListListIOList>("pGrad");

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


TEST_F(SampledPGradTest, CheckInterpolationWorks)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    createWallModelSubregistry(mesh, patch);
    createSamplingHeightField(mesh);

    createPGradField(mesh);
    volVectorField & pGrad = mesh.lookupObjectRef<volVectorField>("pGrad");
    pGrad.primitiveFieldRef() = mesh.C();

    SampledPGradField sampledField(patch, "pointMVC");

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
            { // Here we just check that pointMVC gives us a different value
              // Than stored in the pGrad field cell centres
                ASSERT_NE
                (
                    sampledValues[i][j],
                    pGrad[indexList[i]][j]
                );
            }
        }
    }

}