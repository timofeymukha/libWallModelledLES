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

#include "SampledField.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
namespace Foam
{
    defineTypeNameAndDebug(SampledField, 0);
    // defineRunTimeSelectionTable(Sampler, SamplerRTSTable);
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SampledField::interpolate
(
    scalarListList & values,
    const labelList & cells,
    const pointField & points
) const
{
    FatalErrorInFunction
        << "The sampled field " << name() << " does not support "
        << "interpolation at arbitrary cells and points."
        << exit(FatalError);
}


void Foam::SampledField::postProcess(scalarListList & values) const
{
    FatalErrorInFunction
        << "The sampled field " << name() << " does not support "
        << "post-processing of interpolated values."
        << exit(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //