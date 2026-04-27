#include "codeRules.H"
#include "fvCFD.H"
#include "EquilibriumODEExplicitLawOfTheWall.H"
#include "SingleCellSampler.H"
#include "TOMS748RootFinder.H"
#include "scalarListIOList.H"
#undef Log
#include "gtest.h"
#include "fixtures.H"
#include <functional>

namespace
{
    scalar equilibriumODEUPlus
    (
        const scalar yPlus,
        const scalar kappa,
        const scalar Aplus
    )
    {
        const label nIntervals = 20000;
        const scalar h = yPlus/nIntervals;

        auto integrand =
            [kappa, Aplus](const scalar yp)
            {
                const scalar nutPlus =
                    kappa*yp*sqr(1 - Foam::exp(-yp/Aplus));

                return 1/(1 + nutPlus);
            };

        scalar sum = integrand(0) + integrand(yPlus);

        for (label i = 1; i < nIntervals; i++)
        {
            const scalar yp = i*h;
            sum += (i % 2 ? 4 : 2)*integrand(yp);
        }

        return sum*h/3;
    }

}


class EquilibriumODEExplicitLawOfTheWallTest : public ::testing::Test
{};

class EquilibriumODEExplicitLawOfTheWallChannelFlowTest : public ChannelFlow
{};


TEST_F(EquilibriumODEExplicitLawOfTheWallTest, ConstructFromConstantsHighRe)
{
    EquilibriumODEExplicitLawOfTheWall law(0.387, 15.2516);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.387);
    ASSERT_DOUBLE_EQ(law.Aplus(), 15.2516);
    ASSERT_EQ(law.approximant(), "highRe");
    ASSERT_DOUBLE_EQ(law.p(), 1.1907235645597454);
    ASSERT_DOUBLE_EQ(law.s(), 218.72498935725267);
}


TEST_F(EquilibriumODEExplicitLawOfTheWallTest, DefaultConstantsSelectClassical)
{
    dictionary dict;
    EquilibriumODEExplicitLawOfTheWall law(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.41);
    ASSERT_DOUBLE_EQ(law.Aplus(), 17);
    ASSERT_EQ(law.approximant(), "classical");
    ASSERT_DOUBLE_EQ(law.p(), 1.2029211749614945);
    ASSERT_DOUBLE_EQ(law.s(), 247.0697345675632);
}


TEST_F(EquilibriumODEExplicitLawOfTheWallTest, OtherConstantsSelectGlobal)
{
    dictionary dict;
    dict.add("kappa", 0.395);
    dict.add("Aplus", 16.0);

    EquilibriumODEExplicitLawOfTheWall law(dict);

    ASSERT_EQ(law.approximant(), "global");
    ASSERT_NEAR
    (
        law.p(),
        0.22127889108172996*0.395 + 0.0052185898765612845*16.0
      + 1.0616974939891426,
        1e-12
    );
    ASSERT_NEAR
    (
        law.s(),
        -130.48555166521916*0.395 + 7.964453719974989*16.0
      + 10.6631049793274,
        1e-12
    );
}


TEST_F
(
    EquilibriumODEExplicitLawOfTheWallChannelFlowTest,
    VariantsMatchImplicitEquilibriumODEWithinPaperTolerance
)
{
    extern argList * mainArgs;
    const argList & args = *mainArgs;
    Time runTime(Foam::Time::controlDictName, args);

    autoPtr<fvMesh> meshPtr = createMesh(runTime);
    const fvMesh & mesh = meshPtr();
    createSamplingHeightField(mesh);
    createVelocityField(mesh);

    const fvPatch & patch = mesh.boundary()["bottomWall"];
    SingleCellSampler sampler
    (
        "SingleCellSampler",
        patch,
        3.0,
        "cell",
        "Crawling",
        "CubeRootVol",
        false
    );
    sampler.sample();

    const scalarListIOList& U =
        sampler.db().lookupObject<scalarListIOList>("U");

    struct TestCase
    {
        const char* approximant;
        scalar kappa;
        scalar Aplus;
        scalar tolerance;
    };

    const TestCase testCases[] =
    {
        {"classical", 0.41, 17, 4e-4},
        {"highRe", 0.387, 15.2516, 4e-4},
        {"global", 0.395, 16, 1e-2}
    };

    const label faceI = 5;
    const scalar nu = 8e-6;
    const scalar u = mag(vector(U[faceI][0], U[faceI][1], U[faceI][2]));
    const scalar y = sampler.h()[faceI];

    for (const TestCase& testCase : testCases)
    {
        dictionary dict;
        dict.add("kappa", testCase.kappa);
        dict.add("Aplus", testCase.Aplus);
        dict.add("approximant", word(testCase.approximant));

        EquilibriumODEExplicitLawOfTheWall explicitLaw(dict);

        ASSERT_EQ(explicitLaw.approximant(), word(testCase.approximant));

        const scalar explicitUTau = explicitLaw.uTau(sampler, faceI, nu);

        std::function<scalar(scalar)> value =
            [u, y, nu, &testCase](const scalar uTau)
            {
                return u/uTau
                  - equilibriumODEUPlus
                    (
                        y*uTau/nu,
                        testCase.kappa,
                        testCase.Aplus
                    );
            };

        std::function<scalar(scalar)> derivative =
            [](const scalar)
            {
                return 0;
            };

        TOMS748RootFinder rootFinder("TOMS748", value, derivative, 100);
        const scalar implicitUTau =
            rootFinder.root(explicitUTau, 0.001, 20.0).first;

        const scalar relativeError =
            mag(explicitUTau - implicitUTau)/implicitUTau;

        ASSERT_LT(relativeError, testCase.tolerance);
    }
}
