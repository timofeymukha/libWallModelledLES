#include "fixtures.H"

void createSamplingHeightField(const Foam::fvMesh & mesh)
{
    mesh.time().store
    (
        new volScalarField
        (
            IOobject
            (
                "h",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

Foam::autoPtr<Foam::fvMesh> createMesh(const Foam::Time & runTime)
{
    Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

    meshPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );

    return meshPtr;
}

