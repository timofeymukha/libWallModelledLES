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
    Foam::CellIndexList

Description
    Class for searching for the cells that will be used for sampling.

Authors
    Timofey Mukha, Saleh Rezaeiravesh

SourceFiles
    cellIndexList.C

\*---------------------------------------------------------------------------*/

#ifndef cellIndexList_H
#define cellIndexList_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CellIndexList Declaration
\*---------------------------------------------------------------------------*/

class CellIndexList
{
protected:
    
    // Protected data
        
        //- The patch to build the list for
        const fvPatch & patch_;    
    
        //- The indices of the cells that data is sampled from
        labelList indexList_;
        
        //- A list of wall-normal length-scales associated with the cells
        scalarField lengthList_;
        
        //- The distance from the wall that data is sampled from
        scalarField h_;
                
    // Protected Member Functions
        
        //- Create h field and add to registry
        void createFields() const;
             
        //- Create list of cell-indices from where data is sampled
        void createIndexList();
        
        //- Compute the length-scales
        void createLengthList();
        
public:


    // Constructors

        //- Construct from patch
        CellIndexList
        (
            const fvPatch&
        );

    // Member functions
      
        //- Return the list of cell-indices that are used to sample data
        labelList indexList() const
        {
            return indexList_;
        }
        
        //- Return h
        scalarField h() const
        {
            return h_;
        }
                
        //- Return the length-list
        scalarField lengthList() const
        {
            return lengthList_;
        }
        
        //- Return the patch
        const fvPatch & patch() const
        {
            return patch_;
        }
        
        //- Element access operator
        inline label operator[](const label i) const
        {
            return indexList_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

