#include "codeRules.H"
#include "fvCFD.H"
#include "NewtonRootFinder.H"
#include "SingleCellSampler.H"
#include "SpaldingExplicitLawOfTheWall.H"
#include "SpaldingLawOfTheWall.H"
#undef Log
#include "gtest.h"
#include "fixtures.H"

class SpaldingExplicitLawOfTheWallTest : public ::testing::Test
{};

class SpaldingExplicitLawOfTheWallChannelFlowTest : public ChannelFlow
{};

TEST_F(SpaldingExplicitLawOfTheWallTest, ConstructFromConstantsHighRe)
{
    SpaldingExplicitLawOfTheWall law(0.387, 4.21);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.387);
    ASSERT_DOUBLE_EQ(law.B(), 4.21);
    ASSERT_EQ(law.approximant(), "highRe");
    ASSERT_DOUBLE_EQ(law.p(), 1.1768259148434161);
    ASSERT_DOUBLE_EQ(law.s(), 200.00359862759768);
}


TEST_F(SpaldingExplicitLawOfTheWallTest, DefaultConstantsSelectClassical)
{
    dictionary dict;
    SpaldingExplicitLawOfTheWall law(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.4);
    ASSERT_DOUBLE_EQ(law.B(), 5.5);
    ASSERT_EQ(law.approximant(), "classical");
    ASSERT_DOUBLE_EQ(law.p(), 1.1419065322343513);
    ASSERT_DOUBLE_EQ(law.s(), 400.0);
}


TEST_F(SpaldingExplicitLawOfTheWallTest, OtherConstantsSelectGlobal)
{
    dictionary dict;
    dict.add("kappa", 0.395);
    dict.add("B", 4.0);

    SpaldingExplicitLawOfTheWall law(dict);

    ASSERT_EQ(law.approximant(), "global");
    ASSERT_NEAR(law.p(), 0.1161*0.395 + 0.0203*4.0 + 1.1276, 1e-12);
    ASSERT_NEAR(law.s(), -994.2937*0.395 + 80.3579*4.0 + 301.0408, 1e-12);
}


TEST_F(SpaldingExplicitLawOfTheWallTest, ExplicitApproximantOverridesAuto)
{
    dictionary dict;
    dict.add("kappa", 0.4);
    dict.add("B", 5.5);
    dict.add("approximant", "global");

    SpaldingExplicitLawOfTheWall law(dict);

    ASSERT_EQ(law.approximant(), "global");
    ASSERT_NEAR(law.p(), 0.1161*0.4 + 0.0203*5.5 + 1.1276, 1e-12);
    ASSERT_NEAR(law.s(), -994.2937*0.4 + 80.3579*5.5 + 301.0408, 1e-12);
}


TEST_F
(
    SpaldingExplicitLawOfTheWallChannelFlowTest,
    VariantsMatchImplicitSpaldingWithinPaperTolerance
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

    struct TestCase
    {
        const char* approximant;
        scalar kappa;
        scalar B;
        scalar tolerance;
    };

    const TestCase testCases[] =
    {
        // Fixed three-Gaussian approximants have <0.04% maximum relative error.
        {"classical", 0.4, 5.5, 4e-4},
        {"highRe", 0.387, 4.21, 4e-4},

        // The global one-Gaussian regression is constructed for <1% error.
        {"global", 0.395, 4.8, 1e-2}
    };

    const label faceI = 5;
    const scalar nu = 8e-6;

    for (const TestCase& testCase : testCases)
    {
        SpaldingLawOfTheWall implicitLaw(testCase.kappa, testCase.B);

        dictionary dict;
        dict.add("kappa", testCase.kappa);
        dict.add("B", testCase.B);
        dict.add("approximant", word(testCase.approximant));

        SpaldingExplicitLawOfTheWall explicitLaw(dict);

        ASSERT_EQ(explicitLaw.approximant(), word(testCase.approximant));

        const scalar explicitUTau = explicitLaw.uTau(sampler, faceI, nu);

        std::function<scalar(scalar)> value =
            [&implicitLaw, &sampler, faceI, nu](const scalar uTau)
            {
                return implicitLaw.value(sampler, faceI, uTau, nu);
            };

        std::function<scalar(scalar)> derivative =
            [&implicitLaw, &sampler, faceI, nu](const scalar uTau)
            {
                return implicitLaw.derivative(sampler, faceI, uTau, nu);
            };

        NewtonRootFinder rootFinder("Newton", value, derivative, 100);
        const scalar implicitUTau =
            rootFinder.root(explicitUTau, 0.01, 20.0).first;

        const scalar relativeError =
            mag(explicitUTau - implicitUTau)/implicitUTau;

        ASSERT_LT(relativeError, testCase.tolerance);
    }
}
