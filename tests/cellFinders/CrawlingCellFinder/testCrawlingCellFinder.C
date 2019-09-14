#include "fvCFD.H"
#include "CrawlingCellFinder.H"
#include "gtest.h"
#undef Log
#include "gmock/gmock.h"
#include "fixtures.H"


class CrawlingCellFinderTest : public ChannelFlow
{};

TEST_F(CrawlingCellFinderTest, FullConstructor)
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
    h.boundaryFieldRef()[patch.index()] == 2;
    CrawlingCellFinder finder(patch);

    labelList indexList(patch.size(), 0);
    finder.findCellIndices
    (
        indexList,
        h.boundaryFieldRef()[patch.index()]
    );
    Info << indexList; 

}
