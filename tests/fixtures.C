#include "fixtures.H"

void createSamplingHeightField(const Foam::fvMesh & mesh)
{
    mesh.time().store
    (
        new volScalarField
        (
            IOobject
            (
                "hSampler",
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


void createPGradField(const Foam::fvMesh & mesh)
{
    // Grab h to copy bcs from it.
    const volScalarField & h = mesh.lookupObject<volScalarField>("hSampler");
    
    if (!mesh.foundObject<volVectorField>("pGrad"))
    {
        mesh.thisDb().store
        (     
            new volVectorField
            (
                IOobject
                (
                    "pGrad",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "pGrad",
                    dimLength/sqr(dimTime),
                    pTraits<vector>::zero
                ),
                h.boundaryField().types()
            )
        );
    }
}


void createNutField(const Foam::fvMesh & mesh)
{
    mesh.time().store
    (
        new volScalarField
        (
            IOobject
            (
                "nut",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}


void createNuField(const Foam::fvMesh & mesh, const Foam::volScalarField & nut)
{
    mesh.time().store
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("nu", dimLength*dimLength/dimTime, 0.0),
            nut.boundaryField().types()
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
