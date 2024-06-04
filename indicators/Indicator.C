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

#include "Indicator.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "codeRules.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(Indicator, 0);
    defineRunTimeSelectionTable(Indicator, Patch);
}
#endif

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::Indicator> Foam::Indicator::New
(
    const word & IndicatorName,
    const fvPatch & p
)
{
    auto cstrIter =
    PatchConstructorTablePtr_->find(IndicatorName);

    if (cstrIter == PatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Indicator::New(const word&, const fvPatch & p"

        )   << "Unknown Indicator type "
            << IndicatorName << nl << nl
            << "Valid Indicator types are :" << nl
            << PatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(IndicatorName, p);
}

Foam::autoPtr<Foam::Indicator> Foam::Indicator::New
(
    const dictionary & dict,
    const fvPatch & p
)
{
    word IndicatorName =
        dict.lookupOrDefault<word>("type", "SingleCellIndicator");

    return Foam::Indicator::New(IndicatorName, p);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Indicator::createFields() const
{
    if (!mesh_.thisDb().foundObject<volScalarField>("wmlesIndicator"))
    {
        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "wmlesIndicator",
                    mesh_.thisDb().time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 1.0)
            )
        );
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Indicator::Indicator
(
    const fvPatch & p
)
:
    patch_(p),
    mesh_(patch_.boundaryMesh().mesh())
{
    if (debug)
    {
        Info << "Indicator: Constructing from patch" << nl;
    }

    createFields();
}

Foam::Indicator::Indicator
(
    const word & IndicatorName,
    const fvPatch & p
)
:
    Indicator(p)
{
    if (debug)
    {
        Info << "Indicator: Constructing from name and patch" << nl;
    }
}

Foam::Indicator::Indicator(const Indicator & copy)
:
    patch_(copy.patch_),
    mesh_(copy.mesh_)
{
    if (debug)
    {
        Info << "Indicator: Running copy constructor" << nl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Indicator::~Indicator()
{
}

/*
void sigma_indicator(const volTensorField & gradU)
{
//    auto & U = mesh().thisDb().lookupObject<volVectorField>("U");

    // Limiter
    const dimensionedScalar eps0(dimless, SMALL);
    const dimensionedScalar eps2(dimless/sqr(dimTime), SMALL);
    const dimensionedScalar eps4(dimless/pow4(dimTime), SMALL);
    const dimensionedScalar max2(dimless/sqr(dimTime), GREAT);
    const dimensionedTensor maxTen2
    (
        dimless/sqr(dimTime),
        tensor::max
    );
    const dimensionedTensor minTen2
    (
        dimless/sqr(dimTime),
        tensor::min
    );

 //   volTensorField gradU(fvc::grad(U));
    const volTensorField gradUT(gradU.T());
    const volTensorField G(max(min(gradUT & gradU, maxTen2), minTen2));

    // Tensor invariants
    const volScalarField I1(tr(G));
    const volScalarField I2(0.5*(sqr(I1) - tr(G & G)));
    tmp<volScalarField> tI3 = det(G);

    const volScalarField alpha1(max(sqr(I1)/9.0 - I2/3.0, eps4));

    tmp<volScalarField> talpha2 =
        pow3(min(I1, max2))/27.0 - I1*I2/6.0 + 0.5*tI3;

    const volScalarField alpha3
    (
        1.0/3.0
       *acos
        (
            max
            (
                scalar(-1) + eps0,
                min(scalar(1) - eps0, talpha2/pow(alpha1, 3.0/2.0))
            )
        )
    );

    const scalar piBy3 = constant::mathematical::pi/3.0;
    const volScalarField sigma1
    (
        sqrt(max(I1/3.0 + 2.0*sqrt(alpha1)*cos(alpha3), eps2))
    );
    const volScalarField sigma2
    (
        sqrt(max(I1/3.0 - 2.0*sqrt(alpha1)*cos(piBy3 + alpha3), eps2))
    );
    const volScalarField sigma3
    (
        sqrt(max(I1/3.0 - 2.0*sqrt(alpha1)*cos(piBy3 - alpha3), eps2))
    );

    volScalarField D(
       sigma3
       *(sigma1 - sigma2)
       *(sigma2 - sigma3)
       /max(sqr(sigma1), eps2));


}
*/

void Foam::Indicator::compute
(
    const volScalarField & nu,
    const labelList & inds
) const
{

    const auto & cellCells = mesh_.cellCells();

    const auto & nut = mesh().thisDb().lookupObject<volScalarField>("nut");
    auto & indicator =
        const_cast<volScalarField &>
        (
            mesh().thisDb().lookupObject<volScalarField>("wmlesIndicator")
        );

    scalarField & boundaryValues = indicator.boundaryFieldRef()[patch_.index()];
    scalar C1 = 75;
    scalar C2 = 6;

    forAll(indicator, i)
    {
        indicator[i] =  pow(Foam::tanh(C1 * nut[i] / nu[i]), C2);
    }
    scalarField avrg(inds.size(), 0.0);

    forAll(inds, i)
    {
        label ind = inds[i];
        const labelList& neighbour = cellCells[ind];

        avrg[i] = 0;

        forAll(neighbour, ni)
        {
            avrg[i] += indicator[neighbour[ni]];
        }
        avrg[i] /= neighbour.size();
    }

    forAll(inds, i)
    {
        label ind = inds[i];
        indicator[ind] = avrg[i];
    }

    const scalar threshold = 0.1;
    forAll(boundaryValues, i)
    {
        if (indicator[inds[i]] < threshold)
        {
            boundaryValues[i] = 0.0;
        }
        else
        {
            boundaryValues[i] = 1.0;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
