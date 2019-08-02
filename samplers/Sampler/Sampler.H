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
    2016-2018 Timofey Mukha
    2017      Saleh Rezaeiravesh


SourceFiles
    Sampler.C

\*---------------------------------------------------------------------------*/

#ifndef Sampler_H
#define Sampler_H

#include "tmp.H"
#include "fixedValueFvPatchFields.H"
#include "SampledField.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
        
//- Convert sampled listList to a field
template<class Type>
void listListToField(const scalarListList & list, Field<Type> & field)
{
    forAll(list, i)
    {
        Type element;
        forAll(list[i], j)
        {
            element[j] = list[i][j];
        }
        field[i] = element;
    }
}

    
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

    // Protected Member Functions

        //- Create list of cell-indices from where data is sampled
        virtual void createIndexList() = 0;

        //- Compute the length-scales
        virtual void createLengthList(){}

        //- Project a vector field
        void project(vectorField & field) const;

        //- Create fields
        virtual void createFields();

        //- Compute cell indices suitable for searching sampling cells
        tmp<labelField> findSearchCellLabels() const;

        //- Compute distance field
        tmp<volScalarField> distanceField() const;
        
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
            scalar averagingTime
        );

        //- Construct from type, patch and averaging time
        Sampler
        (
            const word & samplerName,
            const fvPatch &,
            scalar averagingTime
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
            scalar averagingTime
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

        //- Recompute fields to be sampled
        void recomputeFields() const;
        
        //- Sample the fields
        virtual void sample() const = 0; 

        //- Add field for sampling
        virtual void addField(SampledField *);

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
    // RTS tables
        
        // RTS table "PatchAndAveragingTime"
        declareRunTimeSelectionTable
        (
            autoPtr,
            Sampler,
            PatchAndAveragingTime,
            (
                const word& samplerName,
                const fvPatch& p,
                scalar averagingTime
            ),
            (samplerName, p, averagingTime)
        );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

