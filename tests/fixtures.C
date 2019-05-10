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


void createVelocityField(const Foam::fvMesh & mesh)
{
    mesh.time().store
    (
        new volVectorField
        (
            IOobject
            (
                "U",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

void createVelocityGradientField(const Foam::fvMesh & mesh)
{
    mesh.time().store
    (
        new volVectorField
        (
            IOobject
            (
                "wallGradU",
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


void createWallModelSubregistry(const fvMesh & mesh, const fvPatch & patch)
{
    if (!mesh.foundObject<objectRegistry>("wallModelSampling"))
    {
        objectRegistry * subObr = new objectRegistry
        (
            IOobject
            (
                "wallModelSampling",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        subObr->store();
    }

    objectRegistry * subObr = new objectRegistry
    (
        IOobject
        (
            patch.name(),
            mesh.time().constant(),
            mesh.subRegistry("wallModelSampling"),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    subObr->store();
}
