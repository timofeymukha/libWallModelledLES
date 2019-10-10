/*---------------------------------------------------------------------------* \
License
    This file is part of libWallModelledLES.

    libWallModelledLES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libWallModelledLES is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with libWallModelledLES. 
    If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CellFinder.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "codeRules.H"
#include "patchDistMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(CellFinder, 0);
    defineRunTimeSelectionTable(CellFinder, Patch);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //  

Foam::autoPtr<Foam::CellFinder> Foam::CellFinder::New 
(
    const word & CellFinderName,
    const fvPatch & p
)
{
    PatchConstructorTable::iterator cstrIter =
    PatchConstructorTablePtr_->find(CellFinderName);

    if (cstrIter == PatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "CellFinder::New(const word&, const fvPatch & p"
            
        )   << "Unknown CellFinder type "
            << CellFinderName << nl << nl
            << "Valid CellFinder types are :" << nl
            << PatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(CellFinderName, p);
}  

Foam::autoPtr<Foam::CellFinder> Foam::CellFinder::New 
(
    const dictionary & dict,
    const fvPatch & p
)
{
    word CellFinderName =
        dict.lookupOrDefault<word>("type", "SingleCellCellFinder");

    return Foam::CellFinder::New(CellFinderName, p);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::CellFinder::createFields() const
{

}


Foam::tmp<Foam::volScalarField> Foam::CellFinder::distanceField() const
{
    
    // Grab h for the current patch
    const volScalarField & h = mesh_.lookupObject<volScalarField> ("h");
    if (debug)
    {
        Info<< "CellFinder: Creating dist field" << nl;
    }

    tmp<volScalarField> dist
    ( 
        new volScalarField
        (
            IOobject
            (
                "dist",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("dist", dimLength, 0),
            h.boundaryField().types()
        )
    );

    bool precomputedDist = 
    #ifdef FOAM_NEW_GEOMFIELD_RULES
        mag(max(dist().primitiveField())) > VSMALL;
    #else
        mag(max(dist().internalField())) > VSMALL;
    #endif
    
    if (debug)
    {
        Info<<"CellFinder: using precumputed distance field" << nl;
    }

    if (!precomputedDist)
    {
        labelHashSet patchIDs(1);
        patchIDs.insert(patch().index());

        dictionary methodDict = dictionary();
        methodDict.lookupOrAddDefault(word("method"), word("meshWave"));

        if (debug)
        {
            Info<< "Initializing patchDistanceMethod" << nl;
        }

        autoPtr<patchDistMethod> pdm
        (
            patchDistMethod::New
            (
                methodDict,
                mesh_,
                patchIDs
            )
        );

        if (debug)
        {
            Info<< "CellFinder: Computing dist field" << nl;
        }
#ifdef FOAM_NEW_TMP_RULES
        pdm->correct(dist.ref());
#else
        pdm->correct(dist());
#endif
    }

    return dist;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CellFinder::CellFinder
(
    const fvPatch & p
)
:
    patch_(p),
    mesh_(patch_.boundaryMesh().mesh())
{
    if (debug)
    {
        Info << "CellFinder: Constructing from patch" << nl;
    }

    createFields();
}

Foam::CellFinder::CellFinder
(
    const word & CellFinderName,
    const fvPatch & p
)
:
    CellFinder(p)
{
    if (debug)
    {
        Info << "CellFinder: Constructing from name and patch" << nl;
    }
}

Foam::CellFinder::CellFinder(const CellFinder & copy)
:
    patch_(copy.patch_),
    mesh_(copy.mesh_)
{
    if (debug)
    {
        Info << "CellFinder: Running copy constructor" << nl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::CellFinder::~CellFinder()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
