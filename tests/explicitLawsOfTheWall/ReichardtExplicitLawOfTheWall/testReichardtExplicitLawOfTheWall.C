#include "codeRules.H"
#include "fvCFD.H"
#include "NewtonRootFinder.H"
#include "ReichardtExplicitLawOfTheWall.H"
#include "ReichardtLawOfTheWall.H"
#include "SingleCellSampler.H"
#undef Log
#include "gtest.h"
#include "fixtures.H"

class ReichardtExplicitLawOfTheWallTest : public ::testing::Test
{};

class ReichardtExplicitLawOfTheWallChannelFlowTest : public ChannelFlow
{};


TEST_F(ReichardtExplicitLawOfTheWallTest, ConstructFromConstantsHighRe)
{
    ReichardtExplicitLawOfTheWall law(0.387, 11, 3, 6.663050609696008);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.387);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 6.663);
    ASSERT_EQ(law.approximant(), "highRe");
    ASSERT_DOUBLE_EQ(law.p(), 1.2490442812025897);
    ASSERT_DOUBLE_EQ(law.s(), 131.58255333367103);
}


TEST_F(ReichardtExplicitLawOfTheWallTest, DefaultConstantsSelectClassical)
{
    dictionary dict;
    ReichardtExplicitLawOfTheWall law(dict);

    ASSERT_DOUBLE_EQ(law.kappa(), 0.41);
    ASSERT_DOUBLE_EQ(law.B1(), 11);
    ASSERT_DOUBLE_EQ(law.B2(), 3);
    ASSERT_DOUBLE_EQ(law.C(), 7.8);
    ASSERT_EQ(law.approximant(), "classical");
    ASSERT_DOUBLE_EQ(law.p(), 1.2331111626598046);
    ASSERT_DOUBLE_EQ(law.s(), 130.88152902954403);
}


TEST_F(ReichardtExplicitLawOfTheWallTest, OtherConstantsSelectGlobal)
{
    dictionary dict;
    dict.add("kappa", 0.395);
    dict.add("B1", 11);
    dict.add("B2", 3);
    dict.add("C", 7.2);

    ReichardtExplicitLawOfTheWall law(dict);

    ASSERT_EQ(law.approximant(), "global");
    ASSERT_NEAR
    (
        law.p(),
        0.3472185702403766*0.395 - 0.004391948001728183*7.2
      + 1.1436513389571419,
        1e-12
    );
    ASSERT_NEAR
    (
        law.s(),
        -99.55156119159305*0.395 + 19.647991075591285*7.2
      + 18.173005459510563,
        1e-12
    );
}


TEST_F
(
    ReichardtExplicitLawOfTheWallChannelFlowTest,
    VariantsMatchImplicitReichardtWithinPaperTolerance
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
        scalar C;
        scalar tolerance;
    };

    const TestCase testCases[] =
    {
        {"classical", 0.41, 7.8, 4e-4},
        {"highRe", 0.387, 6.663, 4e-4},
        {"global", 0.395, 7.2, 1e-2}
    };

    const label faceI = 5;
    const scalar nu = 8e-6;

    for (const TestCase& testCase : testCases)
    {
        ReichardtLawOfTheWall implicitLaw(testCase.kappa, 11, 3, testCase.C);

        dictionary dict;
        dict.add("kappa", testCase.kappa);
        dict.add("B1", 11);
        dict.add("B2", 3);
        dict.add("C", testCase.C);
        dict.add("approximant", word(testCase.approximant));

        ReichardtExplicitLawOfTheWall explicitLaw(dict);

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
            rootFinder.root(explicitUTau, 0.001, 20.0).first;

        const scalar relativeError =
            mag(explicitUTau - implicitUTau)/implicitUTau;

        ASSERT_LT(relativeError, testCase.tolerance);
    }
}
