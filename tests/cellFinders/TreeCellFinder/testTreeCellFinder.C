#include "fvCFD.H"
#include "TreeCellFinder.H"
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"
#include "fixtures.H"


class TreeCellFinderTest : public ChannelFlow
{};

TEST_F(TreeCellFinderTest, FullConstructor)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);

    volScalarField & h = const_cast<volScalarField &>
    (
        mesh.thisDb().lookupObject<volScalarField>("h")
    );
    const fvPatch & patch = mesh.boundary()["bottomWall"];
    h.boundaryFieldRef()[patch.index()] == 0.25;
    TreeCellFinder finder(patch);

    labelList indexList(patch.size(), 0);
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );
    Info << indexList; 

}
