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

Class
    Foam::SingleCellSampler

@brief
    

Contributors/Copyright:
    2019 Timofey Mukha

SourceFiles
    scalarListListIOList.C

\*---------------------------------------------------------------------------*/

#ifndef scalarListListIOList_H
#define scalarListListIOList_H

#include "scalarList.H"
#include "IOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef List<scalarListList> scalarListListList;
    typedef IOList<scalarListList> scalarListListIOList;
}

#endif
