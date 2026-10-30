#include "dictionary.H"
#include "IFstream.H"
#include "scalarField.H"
#undef Log
#include "gtest.h"
#include <cstdlib>
#include <iostream>
#include <string>

// Runs the channel flow in serial and in parallel, where the sampling cells are
// located on another processor than the wall faces, and checks that the
// results are the same.

// The number of cells and bottomWall faces in the channel flow mesh
const Foam::label nCells = 90;
const Foam::label nWallFaces = 9;

class ParallelSamplerChannelFlow : public ::testing::Test
{
    public:
        ParallelSamplerChannelFlow()
        {
            auto code = std::system("cp -r ../testCases/channel_flow/system .");
            code = std::system("cp -r ../testCases/channel_flow/constant .");
            code = std::system("cp -r ../testCases/channel_flow/0 .");
        }

        virtual ~ParallelSamplerChannelFlow()
        {
            auto code = std::system("rm -r system");
            code = std::system("rm -r constant");
            code = std::system("rm -r 0");
            code = std::system("rm -r 0.01");
            code = std::system("rm -r serial");
            code = std::system("rm -r processor*");
        }
};

class ParallelSamplerIntegrationTest : public ParallelSamplerChannelFlow
{};


// Run a shell command and check that it succeeded
static void run(const std::string & command)
{
    int success = std::system(command.c_str());
    ASSERT_EQ(WIFEXITED(success), true) << command;
    ASSERT_EQ(WEXITSTATUS(success), 0) << command;
}


// Read the internal field, or the value at a patch, of a written field
static Foam::scalarField readValues
(
    const std::string & file,
    const Foam::word & patchName,
    const Foam::label size
)
{
    Foam::IFstream is(file);
    const Foam::dictionary dict(is);

    if (patchName.empty())
    {
        return Foam::scalarField("internalField", dict, size);
    }

    return Foam::scalarField
    (
        "value",
        dict.subDict("boundaryField").subDict(patchName),
        size
    );
}


static void compareValues
(
    const std::string & field,
    const Foam::word & patchName,
    const Foam::label size
)
{
    const Foam::scalarField serial
    (
        readValues("serial/" + field, patchName, size)
    );
    const Foam::scalarField parallel
    (
        readValues("0.01/" + field, patchName, size)
    );

    ASSERT_EQ(serial.size(), size);
    ASSERT_EQ(parallel.size(), size);

    forAll(serial, i)
    {
        ASSERT_NEAR(parallel[i], serial[i], 1e-5*Foam::mag(serial[i]))
            << field << " " << patchName << " " << i;
    }
}


// Run serial and parallel, and compare the sampling cells, the sampling
// distance and the predicted friction velocity. The domain is decomposed in
// the wall-normal direction.
static void runAndCompare
(
    const std::string & setNutDict,
    const int nProcs = 2
)
{
    const std::string n(std::to_string(nProcs));

    ASSERT_NO_FATAL_FAILURE(run("changeDictionary -dict system/" + setNutDict));
    ASSERT_NO_FATAL_FAILURE(run("setExprFields -dict system/setUProfile"));

    ASSERT_NO_FATAL_FAILURE(run("pimpleFoam"));
    ASSERT_NO_FATAL_FAILURE(run("mv 0.01 serial"));

    ASSERT_NO_FATAL_FAILURE
    (
        run
        (
            "foamDictionary -entry numberOfSubdomains -set " + n
          + " system/decomposeParDict"
        )
    );
    ASSERT_NO_FATAL_FAILURE
    (
        run
        (
            "foamDictionary -entry simpleCoeffs/n -set \"(1 " + n + " 1)\""
          + " system/decomposeParDict"
        )
    );
    ASSERT_NO_FATAL_FAILURE(run("decomposePar -force"));
    ASSERT_NO_FATAL_FAILURE(run("mpirun -np " + n + " pimpleFoam -parallel"));
    ASSERT_NO_FATAL_FAILURE(run("reconstructPar"));

    compareValues("samplingCells", "", nCells);
    compareValues("hSampler", "bottomWall", nWallFaces);
    compareValues("uTauPredicted", "bottomWall", nWallFaces);
}


TEST_F(ParallelSamplerIntegrationTest, CrawlingIndexSerialEqualsParallel)
{
    runAndCompare("setNutLOTWSpaldingParallelCrawlingIndex");
}


TEST_F(ParallelSamplerIntegrationTest, CrawlingDistanceSerialEqualsParallel)
{
    runAndCompare("setNutLOTWSpaldingParallelCrawlingDistance");
}


TEST_F(ParallelSamplerIntegrationTest, NormalRaySerialEqualsParallel)
{
    runAndCompare("setNutLOTWSpaldingParallelNormalRay");
}


// With four processors, the crawls are handed over several times, and some
// processors hold neither wall faces nor sampling cells
TEST_F(ParallelSamplerIntegrationTest, CrawlingIndexSerialEqualsFourProcessors)
{
    runAndCompare("setNutLOTWSpaldingParallelCrawlingIndex", 4);
}


TEST_F(ParallelSamplerIntegrationTest, NormalRaySerialEqualsFourProcessors)
{
    runAndCompare("setNutLOTWSpaldingParallelNormalRay", 4);
}
