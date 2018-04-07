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

#include "Sampler.H"
#include "meshSearch.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "objectRegistry.H"
#include "IOField.H"
#include "SampledField.H"
#include "SampledVelocityField.H"
#include "SampledPGradField.H"
#include "SampledWallGradUField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Sampler::createFields()
{
   
    //const objectRegistry & db = patch_.boundaryMesh().mesh();
    
    if (!mesh_.thisDb().found("h"))
    {
        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch_.boundaryMesh().mesh()
            )
        );
    }
    
    volScalarField & h = 
        const_cast<volScalarField &>(mesh_.lookupObject<volScalarField> ("h"));
    
   // Field that marks cells that are used for sampling
    if (!mesh_.thisDb().found("samplingCells"))
    {
        mesh_.thisDb().store
        (
            new volScalarField
            (
                IOobject
                (
                    "samplingCells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("samplingCells", dimless,0),
                h.boundaryField().types()
            )
        );
    }
      
}


void Foam::Sampler::createIndexList()
{
   
    const label patchIndex = patch().index();
    
    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh();
    
    // Grab h for the current patch
    volScalarField & h = 
        const_cast<volScalarField &>(mesh.lookupObject<volScalarField> ("h"));
    
    h_ = h.boundaryField()[patchIndex];



    // Create a searcher for the mesh
    meshSearch ms(mesh);
    
    // Grab face centres, normal and adjacent cells' centres to each patch face
    const vectorField & faceCentres = patch().Cf();
    const tmp<vectorField> tfaceNormals = patch().nf();
    const vectorField faceNormals = tfaceNormals();
    const tmp<vectorField> tcellCentres = patch().Cn();
    const vectorField cellCentres = tcellCentres();
    
    // Grab the global indices of adjacent cells 
    const labelUList & faceCells = patch().faceCells();

    vector point;
    forAll(faceCentres, i)
    {
        // If h is zero, set it to distance to adjacent cell's centre
        // Set the cellIndexList component accordingly.
        if (h_[i] == 0)
        {          
            h_[i] = mag(cellCentres[i] - faceCentres[i]);
            indexList_[i] = faceCells[i];          
        }
        else
        {
            // Grab the point h away along the face normal
            point = faceCentres[i] - faceNormals[i]*h_[i];

            // Check that point is inside the (processor) domain
            // Otherwise fall back to adjacent cell's centre.
            if (!ms.isInside(point))
            {
                point = cellCentres[i];
            }

            // Find the cell where the point is located
            indexList_[i] = ms.findNearestCell(point, -1, true);
                        
            // Set h to the distance between face centre and located cell's
            // center
            h_[i] = mag(mesh.C()[indexList_[i]] - faceCentres[i]);
        }
    }
    
    // Assign computed h_ to the global h field
    h.boundaryField()[patch().index()] == h_;
    
    // Grab samplingCells field
    volScalarField & samplingCells = 
        const_cast<volScalarField &>
        (
            mesh.lookupObject<volScalarField> ("samplingCells")
        );
    
    forAll(indexList_, i)
    {
        samplingCells[indexList_[i]] = patchIndex; 
    }
}


void Foam::Sampler::createLengthList()
{
    // Grab the mesh
    const fvMesh & mesh = patch().boundaryMesh().mesh(); 
    
    // Cell volumes
    const scalarField & V = mesh.V();
    
    forAll(lengthList_, i)
    {
        lengthList_[i] = pow(V[indexList_[i]], 1.0/3.0);
    }
    
}

void Foam::Sampler::project(vectorField & field) const
{
    const tmp<vectorField> tfaceNormals = patch_.nf();
    const vectorField faceNormals = tfaceNormals();

    
    forAll(field, i)
    {    
        // Normal component as dot product with (inwards) face normal
        vector normal = -faceNormals[i]*(field[i] & -faceNormals[i]);
        
        // Subtract normal component to get the parallel one
        field[i] -= normal;        
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Sampler::Sampler
(
    const fvPatch & p,
    scalar averagingTime
)
:
    patch_(p),
    indexList_(p.size()),
    lengthList_(p.size()),
    h_(p.size(), 0),
    averagingTime_(averagingTime),
    db_(patch_.boundaryMesh().mesh().subRegistry("wallModelSampling", 1)),
    mesh_(patch_.boundaryMesh().mesh()),
    sampledFields_(0)
{
    createFields();
    createIndexList();
    createLengthList();
    
    addField
    (
            new SampledVelocityField(patch_, indexList_)     
    );
    
    addField
    (
            new SampledWallGradUField(patch_, indexList_)     
    );
    
    //recomputeFields();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Sampler::~Sampler()
{
    forAll(sampledFields_, i)
    {
        delete sampledFields_[i];
    }
}

void Foam::Sampler::sample() const
{

    // Weight for time-averaging, default to 1 i.e no averaging.
    scalar eps = 1;
    if (averagingTime_ > mesh_.time().deltaTValue())
    {
        eps = mesh_.time().deltaTValue()/averagingTime_;
    }
    
    Info << db_.names() << nl;
    
    forAll(sampledFields_, fieldI)
    {

        List<List<scalar> > sampledList = sampledFields_[fieldI]->sample();
        //Info << "sampledList " << sampledList << nl;
        
        if (sampledFields_[fieldI]->nDims() == 3)
        {
            //const vectorField & sampledField = 
            //    listListToField<vector>(sampledList);
            vectorField sampledField(patch_.size());
            listListToField<vector>(sampledList, sampledField);
            //Info << "sampledField " << sampledField << nl;
            
            vectorField & storedValues = const_cast<vectorField &>
            (
                db_.lookupObject<vectorField>(sampledFields_[fieldI]->name())
            );
            //Info << "storedvalues" << storedValues << nl;
            storedValues = eps*sampledField + (1 - eps)*storedValues;
        }
        

        //sampled[i] = eps*sampleField[i] + (1 - eps)*sampled[i];
    }
    
}

template<class Type>
//Foam::tmp<Foam::Field<Type> >
//Foam::Field<Type>
//Foam::Sampler::listListToField(const List<List<scalar> > & list) const
void Foam::Sampler::listListToField
(
    const List<List<scalar> > & list,
    Field<Type> & field    
) const
{
//    tmp<Field<Type> > tField(new Field<Type>(list.size()));
//    Field<Type> & field = tField();

//    Field<Type> field(list.size());
    scalar nDims = list[0].size();
    
    forAll(list, i)
    {
        Type element;
        for (int j=0; j<nDims; j++)
        {
            element[j] = list[i][j];
        }
        field[i] = element;
        
    }
    //return tField;
    //return field;
}



void Foam::Sampler::addField(SampledField * field)
{
    sampledFields_.append(field);
}


void Foam::Sampler::recomputeFields() const
{
    forAll(sampledFields_, i)
    {
        sampledFields_[i]->recompute();
    } 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //