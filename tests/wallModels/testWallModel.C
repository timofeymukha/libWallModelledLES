#include "codeRules.H"
#include "fvCFD.H"
#include "wallModelFvPatchScalarField.H"
#include "directFvPatchFieldMapper.H"
#undef Log
#include "gtest.h"
#include "gmock/gmock.h"
#include "fixtures.H"


namespace Foam
{
    class DummyWallModel : public wallModelFvPatchScalarField
    {
        public:
            TypeName("DummyWallModel");

            DummyWallModel
            (
                const fvPatch& p,
                const DimensionedField<scalar, volMesh>& iF
            )
            :
            wallModelFvPatchScalarField(p, iF)
            {}

            DummyWallModel
            (
                const fvPatch& p,
                const DimensionedField<scalar, volMesh>& f,
                const dictionary& d
            )
            :
            wallModelFvPatchScalarField(p, f, d)
            {}

            DummyWallModel
            (
                const DummyWallModel& orig,
                const fvPatch& p,
                const DimensionedField<scalar, volMesh>& f,
                const fvPatchFieldMapper& mapper
            )
            :
            wallModelFvPatchScalarField(orig, p, f, mapper)
            {}

            DummyWallModel
            (
                const DummyWallModel& orig,
                const DimensionedField<scalar, volMesh>& f
            )
            :
            wallModelFvPatchScalarField(orig, f)
            {}

            DummyWallModel(const DummyWallModel &) = default;

            
        // Destructor
            virtual ~DummyWallModel()
            {}

            virtual tmp<scalarField> calcNut() const
            {
                return tmp<scalarField>(new scalarField(patch().size(), 0.0));
            }
            
    };

    makePatchTypeField
    (
        fvPatchScalarField,
        DummyWallModel
    );
}

class WallModelTest : public ChannelFlow
{

    public:
        WallModelTest()
        :
        ChannelFlow()
        {
            system("cp 0/nutFixedValue 0/nut");
        }
};


TEST_F(WallModelTest, ConstructorW1)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    const volScalarField nutField = mesh.lookupObject<volScalarField>("nut");

    DummyWallModel model(patch, nutField);

    ASSERT_EQ(&patch, &model.patch());
    ASSERT_EQ(&nutField, &model.internalField());
    ASSERT_FLOAT_EQ(model.averagingTime(), 0.0);
    ASSERT_FLOAT_EQ(model.consumedTime(), 0.0);
    ASSERT_EQ(model.copyToPatchInternalField(), false);
    ASSERT_TRUE(mesh.foundObject<volScalarField>("hSampler"));
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallShearStress"));
    ASSERT_TRUE(mesh.foundObject<volScalarField>("uTauPredicted"));
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallGradU"));
}

TEST_F(WallModelTest, ConstructorW2)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];

    labelList addressing(patch.size());
    
    // we will map to self
    forAll(patch, i)
    {
        addressing[i] = i;
    }

    directFvPatchFieldMapper mapper(addressing);

    const volScalarField nutField = mesh.lookupObject<volScalarField>("nut");

    DummyWallModel model(patch, nutField);
    DummyWallModel model2(model, patch, nutField, mapper);

    ASSERT_EQ(&patch, &model.patch());
    ASSERT_EQ(&nutField, &model.internalField());
    ASSERT_FLOAT_EQ(model2.averagingTime(), 0.0);
    ASSERT_FLOAT_EQ(model2.consumedTime(), 0.0);
    ASSERT_EQ(model2.copyToPatchInternalField(), false);
}

TEST_F(WallModelTest, ConstructorW3)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);

    dictionary dict;
    dict.add("averagingTime", 0.1);
    dict.add("copyToPatchInternalField", true);
    dict.add("value", "uniform 0.0");

    const volScalarField nutField = mesh.lookupObject<volScalarField>("nut");
    const fvPatch & patch = mesh.boundary()["bottomWall"];

    DummyWallModel model(patch, nutField, dict);
    ASSERT_FLOAT_EQ(model.averagingTime(), 0.1);
    ASSERT_FLOAT_EQ(model.consumedTime(), 0.0);
    ASSERT_EQ(model.copyToPatchInternalField(), true);

    ASSERT_TRUE(mesh.foundObject<volScalarField>("hSampler"));
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallShearStress"));
    ASSERT_TRUE(mesh.foundObject<volScalarField>("uTauPredicted"));
    ASSERT_TRUE(mesh.foundObject<volVectorField>("wallGradU"));
}

TEST_F(WallModelTest, CopyConstructorW4)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);

    dictionary dict;
    dict.add("averagingTime", 0.1);
    dict.add("value", "uniform 0.0");
    dict.add("copyToPatchInternalField", true);

    const volScalarField nutField = mesh.lookupObject<volScalarField>("nut");
    const fvPatch & patch = mesh.boundary()["bottomWall"];

    DummyWallModel model(patch, nutField, dict);
    DummyWallModel model2(model);


    ASSERT_DOUBLE_EQ(model2.averagingTime(), 0.1);
    ASSERT_DOUBLE_EQ(model2.consumedTime(), 0.0);
    ASSERT_EQ(model2.copyToPatchInternalField(), true);
}

TEST_F(WallModelTest, CopyConstructorW5)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);

    dictionary dict;
    dict.add("averagingTime", 0.1);
    dict.add("value", "uniform 0.0");

    const volScalarField nutField = mesh.lookupObject<volScalarField>("nut");
    const fvPatch & patch = mesh.boundary()["bottomWall"];

    DummyWallModel model(patch, nutField, dict);
    DummyWallModel model2(model, nutField);

    ASSERT_DOUBLE_EQ(model2.averagingTime(), 0.1);
    ASSERT_DOUBLE_EQ(model2.consumedTime(), 0.0);
    ASSERT_EQ(model2.copyToPatchInternalField(), false);
}

TEST_F(WallModelTest, UpdateCoeffs)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);
    const volScalarField & nutField =  mesh.lookupObject<volScalarField>("nut");
    createNuField(mesh, nutField);
    createVelocityField(mesh);

    dictionary dict;
    dict.add("averagingTime", 0.1);
    dict.add("value", "uniform 0.0");

    const fvPatch & patch = mesh.boundary()["bottomWall"];
        
    DummyWallModel model(patch, nutField, dict);
    volScalarField & nuField = const_cast<volScalarField &>
    (
        mesh.lookupObject<volScalarField>("nu")
    );

    volVectorField & wallGradUField = const_cast<volVectorField &>
    (
        mesh.lookupObject<volVectorField>("wallGradU")
    );

    volVectorField & wssField = const_cast<volVectorField &>
    (
        mesh.lookupObject<volVectorField>("wallShearStress")
    );
    nuField.boundaryFieldRef()[patch.index()] == 2;
    wallGradUField.boundaryFieldRef()[patch.index()] == vector(1, 2, 3);

    model.wallModelFvPatchScalarField::updateCoeffs();
    
    const vectorField & wss = wssField.boundaryFieldRef()[patch.index()];

    forAll(wss, i)
    {
        ASSERT_DOUBLE_EQ(wss[i][0], 2);
        ASSERT_DOUBLE_EQ(wss[i][1], 4);
        ASSERT_DOUBLE_EQ(wss[i][2], 6);
    }
    ASSERT_DOUBLE_EQ(model.averagingTime(), 0.1);
    ASSERT_DOUBLE_EQ(model.copyToPatchInternalField(), false);
}      