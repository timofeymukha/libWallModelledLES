#include "codeRules.H"
#include "fvCFD.H"
#include "wallModelFvPatchScalarField.H"
#include "directFvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
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
                return tmp<scalarField>(new scalarField(patch().size(), 2.0));
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
    ASSERT_EQ(model.silent(), false);
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
    ASSERT_EQ(model.silent(), false);

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
    ASSERT_EQ(model2.silent(), false);
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
    ASSERT_EQ(model2.silent(), false);
}

TEST_F(WallModelTest, ConstructorPreservesRegisteredHSampler)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);
    createSamplingHeightField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    volScalarField & h = mesh.lookupObjectRef<volScalarField>("hSampler");
    h.boundaryFieldRef()[patch.index()] == 42.0;

    const volScalarField & nutField = mesh.lookupObject<volScalarField>("nut");

    DummyWallModel model(patch, nutField);

    const volScalarField & hAfter =
        mesh.lookupObject<volScalarField>("hSampler");

    forAll(hAfter.boundaryField()[patch.index()], i)
    {
        ASSERT_DOUBLE_EQ(hAfter.boundaryField()[patch.index()][i], 42.0);
    }
}

TEST_F(WallModelTest, UpdateCoeffs)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);
    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createNutField(mesh);
    const volScalarField & nutField = mesh.lookupObject<volScalarField>("nut");

    dictionary dict;
    dict.add("averagingTime", 0.1);
    dict.add("copyToPatchInternalField", true);
    dict.add("silent", true);
    dict.add("value", "uniform 0.0");

    const fvPatch & patch = mesh.boundary()["bottomWall"];
        
    DummyWallModel model(patch, nutField, dict);

    model.wallModelFvPatchScalarField::updateCoeffs();

    const scalarField & patchNut = model;
    forAll(patchNut, i)
    {
        ASSERT_DOUBLE_EQ(patchNut[i], 2.0);
    }

    const labelUList faceCells = patch.faceCells();
    forAll(faceCells, i)
    {
        ASSERT_DOUBLE_EQ(nutField[faceCells[i]], 2.0);
    }

    ASSERT_DOUBLE_EQ(model.averagingTime(), 0.1);
    ASSERT_DOUBLE_EQ(model.copyToPatchInternalField(), true);
    ASSERT_EQ(model.silent(), true);
}
