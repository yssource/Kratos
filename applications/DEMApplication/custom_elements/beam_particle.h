//
// Author: Joaquín Irazábal, jirazabal@cimne.upc.edu
//

#if !defined(KRATOS_BEAM_PARTICLE_H_INCLUDED)
#define KRATOS_BEAM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "spheric_continuum_particle.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) BeamParticle : public SphericContinuumParticle {

        public:

        /// Pointer definition of BeamParticle
        KRATOS_CLASS_POINTER_DEFINITION(BeamParticle);

        BeamParticle() : SphericContinuumParticle() {}
        BeamParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry) {}
        BeamParticle(IndexType NewId, NodesArrayType const& ThisNodes) : SphericContinuumParticle(NewId, ThisNodes) {}
        BeamParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SphericContinuumParticle(NewId, pGeometry, pProperties) {}

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return SphericContinuumParticle::Pointer(new BeamParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Destructor
        virtual ~BeamParticle() {};

        /// Turn back information as a string
        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "BeamParticle" ;
            return buffer.str();
        }


        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "BeamParticle";}

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const override {}

        virtual void ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer &,
                                                   ProcessInfo& r_process_info,
                                                   array_1d<double, 3>& rElasticForce,
                                                   array_1d<double, 3>& rContactForce,
                                                   double& RollingResistance) override;

        private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        /*
        /// Assignment operator
        BeamParticle& operator=(BeamParticle const& rOther)
        {
        return *this;
        }

        /// Copy constructor
        BeamParticle(BeamParticle const& rOther)
        {
        *this = rOther;
        }
        */

        ///@}

    }; // Class BeamParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, BeamParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const BeamParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_BEAM_PARTICLE_H_INCLUDED defined
