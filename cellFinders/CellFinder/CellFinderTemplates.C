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
#include "PstreamBuffers.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::CellFinder::exchange
(
    const UList<DynamicList<Type>> & sendLists,
    List<List<Type>> & recvLists
)
{
    const label myProc = Pstream::myProcNo();

    recvLists.setSize(Pstream::nProcs());
    recvLists[myProc] = sendLists[myProc];

    if (!Pstream::parRun())
    {
        return;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(sendLists, proci)
    {
        if (proci != myProc && sendLists[proci].size())
        {
            UOPstream toProc(proci, pBufs);
            toProc << sendLists[proci];
        }
    }

    pBufs.finishedSends();

    forAll(recvLists, proci)
    {
        if (proci != myProc && pBufs.recvDataCount(proci))
        {
            UIPstream fromProc(proci, pBufs);
            fromProc >> recvLists[proci];
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
