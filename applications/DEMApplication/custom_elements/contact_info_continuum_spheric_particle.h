//
// Author: Joaquín Irazábal jirazabal@cimne.upc.edu
//


#if !defined(KRATOS_CONTACT_INFO_CONTINUUM_SPHERIC_PARTICLE_H_INCLUDED)
#define  KRATOS_CONTACT_INFO_CONTINUUM_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "spheric_continuum_particle.h"


namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) ContactInfoContinuumSphericParticle : public SphericContinuumParticle
{
public:

/// Pointer definition of ContactInfoContinuumSphericParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ContactInfoContinuumSphericParticle);

// typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
// typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;

// typedef GlobalPointersVector<Element> ParticleWeakVectorType;
// typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
// typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
typedef SphericContinuumParticle BaseType;
typedef BaseType::ParticleDataBuffer BaseBufferType;
typedef std::unique_ptr<BaseType::ParticleDataBuffer> BaseBufferPointerType;

/// Default constructor.
ContactInfoContinuumSphericParticle();
ContactInfoContinuumSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry);
ContactInfoContinuumSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
ContactInfoContinuumSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
ContactInfoContinuumSphericParticle(Element::Pointer p_continuum_spheric_particle);

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~ContactInfoContinuumSphericParticle(){}

ContactInfoContinuumSphericParticle& operator=(const ContactInfoContinuumSphericParticle& rOther);

/// Turn back information as a string.
std::string Info() const override
{
std::stringstream buffer;
buffer << "ContactInfoContinuumSphericParticle" ;
return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream& rOStream) const override {rOStream << "ContactInfoContinuumSphericParticle";}

/// Print object's data.
void PrintData(std::ostream& rOStream) const override {}

virtual double GetParticleInitialCohesion();
void   SetParticleInitialCohesionFromProperties(double* particle_initial_cohesion);
virtual double GetAmountOfCohesionFromStress();
void   SetAmountOfCohesionFromStressFromProperties(double* amount_of_cohesion_from_stress);
virtual double GetParticleConicalDamageContactRadius();
void   SetParticleConicalDamageContactRadiusFromProperties(double* particle_contact_radius);
virtual double GetParticleConicalDamageMaxStress();
void   SetParticleConicalDamageMaxStressFromProperties(double* particle_max_stress);
virtual double GetParticleConicalDamageGamma();
void   SetParticleConicalDamageGammaFromProperties(double* particle_gamma);
virtual double GetLevelOfFouling();
void   SetLevelOfFoulingFromProperties(double* level_of_fouling);

double SlowGetParticleInitialCohesion();
double SlowGetAmountOfCohesionFromStress();
double SlowGetParticleConicalDamageContactRadius();
double SlowGetParticleConicalDamageMaxStress();
double SlowGetParticleConicalDamageGamma();
double SlowGetLevelOfFouling();

std::vector<double> mNeighbourContactRadius;
std::vector<double> mNeighbourRigidContactRadius;
std::vector<double> mNeighbourIndentation;
std::vector<double> mNeighbourRigidIndentation;
std::vector<double> mNeighbourTgOfFriAng;
std::vector<double> mNeighbourRigidTgOfFriAng;
std::vector<double> mNeighbourContactStress;
std::vector<double> mNeighbourRigidContactStress;
std::vector<double> mNeighbourCohesion;
std::vector<double> mNeighbourRigidCohesion;

protected:

class ParticleDataBuffer: public SphericParticle::ParticleDataBuffer
{
public:

ParticleDataBuffer(SphericParticle* p_this_particle): SphericParticle::ParticleDataBuffer(p_this_particle){}

virtual ~ParticleDataBuffer(){}

};

std::unique_ptr<SphericParticle::ParticleDataBuffer> CreateParticleDataBuffer(SphericParticle* p_this_particle) override
{
    return std::unique_ptr<SphericParticle::ParticleDataBuffer>(new ParticleDataBuffer(p_this_particle));
}

void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                        std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces) override;

void ComputeNewRigidFaceNeighboursHistoricalData() override;

private:

void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
}

void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
}

}; // Class ContactInfoContinuumSphericParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            ContactInfoContinuumSphericParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const ContactInfoContinuumSphericParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_CONTACT_INFO_CONTINUUM_SPHERIC_PARTICLE_H_INCLUDED  defined
