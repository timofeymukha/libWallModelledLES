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
    Class for sampling data to the wall models for a single cell per face.

Contributors/Copyright:
    2018-2021 Timofey Mukha

SourceFiles
    SingleCellSampler.C

\*---------------------------------------------------------------------------*/

#ifndef SingleCellSampler_H
#define SingleCellSampler_H

#include "fixedValueFvPatchFields.H"
#include "Sampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SampledField;

/*---------------------------------------------------------------------------*\
                         Class SingleCellSampler Declaration
\*---------------------------------------------------------------------------*/

class SingleCellSampler: public Sampler
{
protected:
    
    // Protected data
        
        //- The indices of the cells that data is sampled from
        labelList indexList_;
        
        //- A list of wall-normal length-scales associated with the cells
        scalarField lengthList_;
        
        //- The distance from the wall that data is sampled from
        scalarField h_;
        
    // Protected Member Functions
                    
        //- Create list of cell-indices from where data is sampled
        void createIndexList() override;

        //- Compute the length-scales
        void createLengthList(const word lengthScaleType) override;
        
        //- Compute length-scales as cubic root of the volume
        void createLengthListCubeRootVol();
        
        //- Compute length-scales as distance across wall-normal direction
        void createLengthListWallNormalDistance();
        
public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("SingleCellSampler");
#endif

    // Constructors

        SingleCellSampler
        (
            const fvPatch&,
            scalar averagingTime,
            const word interpolationType="cell",
            const word cellFinderType="Tree",
            const word lengthScaleType="CubeRootVol",
            bool hIsIndex=false,
            bool excludeWallAdjacent=false
        );

        SingleCellSampler
        (
            const word & samplerName,
            const fvPatch & p,
            scalar averagingTime,
            const word interpolationType="cell",
            const word cellFinderType="Tree",
            const word lengthScaleType="CubeRootVol",
            bool hIsIndex=false,
            bool excludeWallAdjacent=false
        );
        
        SingleCellSampler(const SingleCellSampler &) = default;
        
    // Destructor
        virtual ~SingleCellSampler();

    // Member functions
        //- Add field for sampling
        void addField(SampledField *) override;

        //- Return the list of cell-indices that are used to sample data
        const labelList & indexList() const
        {
            return indexList_;
        }
        
        //- Return h
        const scalarField & h() const
        {
            return h_;
        }
                
        //- Return the length-list
        virtual const scalarField & lengthList() const
        {
            return lengthList_;
        }
        
        //- Element access operator
        inline label operator[](const label i) const
        {
            return indexList_[i];
        }
        
        //- Sample the fields
        void sample() const override; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

