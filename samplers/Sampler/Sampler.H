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
    Foam::Sampler

Description
    Class for sampling data to the wall models.

Contributors/Copyright:
    2016-2021 Timofey Mukha
    2017      Saleh Rezaeiravesh


SourceFiles
    Sampler.C

\*---------------------------------------------------------------------------*/

#ifndef Sampler_H
#define Sampler_H

#include "tmp.H"
#include "fixedValueFvPatchFields.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SampledField;

/*---------------------------------------------------------------------------*\
                         Class Sampler Declaration
\*---------------------------------------------------------------------------*/

class Sampler
{
protected:
    
    // Protected data
        
        //- The patch to build the list for
        const fvPatch & patch_;    

        //- Time-averaging scale of the sampled values
        scalar averagingTime_;

        //- The global object registry
        const fvMesh & mesh_;

        //- List of sampled fields
        PtrList<SampledField> sampledFields_;

        //- Type of interpolation to do within the sampling cell
        //- cell, cellPoint or cellPointFace. Defaults to cell.
        word interpolationType_;

        //- Type of cellFinder to use for finding sampling cells.
        word cellFinderType_;

        //- Way to compute the length-scale for the cells
        word lengthScaleType_;

        //- Whether the h field is the consecutive index of the samping cell
        bool hIsIndex_;

    // Protected Member Functions

        //- Create list of cell-indices from where data is sampled
        virtual void createIndexList() = 0;

        //- Compute the length-scales
        virtual void createLengthList(const word lengthScaleType){}

        //- Create fields
        virtual void createFields();

public:

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    //- Runtime type information
        TypeName("Sampler");
#endif

    // Constructors

        //- Construct from patch and averaging time
        Sampler
        (
            const fvPatch &,
            scalar averagingTime,
            const word interpolationType,
            const word cellFinderType,
            const word lengthScaleType,
            bool hIsIndex
        );

        //- Construct from type, patch and averaging time
        Sampler
        (
            const word & samplerName,
            const fvPatch &,
            scalar averagingTime,
            const word interpolationType,
            const word cellFinderType,
            const word lengthScaleType,
            bool hIsIndex
        );


        //- Copy constructor
        Sampler(const Sampler &);
        
        //- Destructor
        virtual ~Sampler();

    // Selectors
        static autoPtr<Sampler> New
        (
            const word & samplerName,
            const fvPatch &,
            scalar averagingTime,
            const word interpolationType,
            const word cellFinderType,
            const word lengthScaleType,
            bool hIsIndex
        );

        static autoPtr<Sampler> New
        (
            const dictionary &,
            const fvPatch &
        );

    // Member functions

        //- Return the patch
        const fvPatch & patch() const
        {
            return patch_;
        }
        
        //- Get the registry with sampled fields
        const objectRegistry & db() const
        {
            return mesh().subRegistry("wallModelSampling").subRegistry(patch().name());
        }
        
        //- Get the mesh
        const fvMesh & mesh() const
        {
            return mesh_;
        }

        //- Get averaging time
        scalar averagingTime() const
        {
            return averagingTime_;
        }

        label nSampledFields() const
        {
            return sampledFields_.size();
        }

        word interpolationType() const
        {
            return interpolationType_;
        }

        word cellFinderType() const
        {
            return cellFinderType_;
        }
        
        word lengthScaleType() const
        {
            return lengthScaleType_;
        }

        bool hIsIndex() const
        {
            return hIsIndex_;
        }

        //- Recompute fields to be sampled
        void recomputeFields() const;
        
        //- Sample the fields
        virtual void sample() const = 0; 

        //- Add field for sampling
        virtual void addField(SampledField *);

        //- Write properties to dictionary
        void write(Ostream & os) const;

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    // RTS tables
        
        // RTS table "SamplerRTSTable"
        declareRunTimeSelectionTable
        (
            autoPtr,
            Sampler,
            SamplerRTSTable,
            (
                const word& samplerName,
                const fvPatch& p,
                scalar averagingTime,
                const word interpolationType,
                const word cellFinderType,
                const word lengthScaleType,
                bool hIsIndex
            ),
            (
                samplerName,
                p,
                averagingTime,
                interpolationType,
                cellFinderType,
                lengthScaleType,
                hIsIndex
            )
        );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

