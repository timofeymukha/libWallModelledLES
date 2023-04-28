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
    Foam::MultiCellSampler

@brief
    Class for sampling from several consecutive cells per wall face.

Contributors/Copyright:
    2018-2021 Timofey Mukha

SourceFiles
    MultiCellSampler.C

\*---------------------------------------------------------------------------*/

#ifndef MultiCellSampler_H
#define MultiCellSampler_H

#include "fixedValueFvPatchFields.H"
#include "Sampler.H"
#include "SampledField.H"
#include "scalarListIOList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    
/*---------------------------------------------------------------------------*\
                         Class MultiCellSampler Declaration
\*---------------------------------------------------------------------------*/

class MultiCellSampler: public Sampler
{
protected:
    
    // Protected data
        
        //- The lists of indices of the cells that data is sampled from
        labelListList indexList_;
        
        //- The distance from the wall that data is sampled from
        scalarListList h_;

        //- The distance from the wall that data is sampled from
        scalarListList lengthList_;

        
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
        TypeName("MultiCellSampler");
#endif

    // Constructors

        //- Construct from patch and averaging time
        MultiCellSampler
        (
            const fvPatch&,
            scalar averagingTime,
            const word interpolationType="cell",
            const word cellFinderType="Tree",
            const word lengthScaleType="CubeRootVol",
            bool hIsIndex=false,
            bool excludeWallAdjacent=false
        );

        //- Construct from type, patch and averaging time
        MultiCellSampler
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
        
        //- Copy constructor
        MultiCellSampler(const MultiCellSampler & orig)
        :
        Sampler(orig)
        {}


    // Destructor
        virtual ~MultiCellSampler()
        {}

    // Member functions

        //- Add field for sampling
        virtual void addField(SampledField *);

        //- Return the list of lists of cell-indices that are used to sample data
        const labelListList & indexList() const
        {
            return indexList_;
        }


        //- Return h
        virtual const scalarListList & h() const
        {
            return h_;
        }

                
        //- Return the length-list
        virtual const scalarListList & lengthList() const
        {
            return lengthList_;
        }

        //- Returning sampling cell indices for a given wall face using []
        inline labelList operator[](const label i) const
        {
            return indexList_[i];
        }
        
        //- Sample the fields
        virtual void sample() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

