// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Beam_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Beam::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Beam(*this));
        return p_clone;
    }

    void DEM_KDEM_Beam::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Beam to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM_Beam::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                         SphericContinuumParticle* neighbor,
                                                         double equiv_young,
                                                         double distance,
                                                         double calculation_area,
                                                         double LocalCoordSystem[3][3],
                                                         double ElasticLocalRotationalMoment[3],
                                                         double ViscoLocalRotationalMoment[3],
                                                         double equiv_poisson,
                                                         double indentation) {

        KRATOS_TRY
        double LocalDeltaRotatedAngle[3]    = {0.0};
        double LocalDeltaAngularVelocity[3] = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3> GlobalDeltaAngularVelocity;
        noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);

        const double element_mass  = element->GetMass();
        const double neighbor_mass = neighbor->GetMass();
        const double equiv_mass    = element_mass * neighbor_mass / (element_mass + neighbor_mass);
        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double Inertia_Ix = std::max(element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X], neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X]);
        const double Inertia_Iy = std::max(element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y], neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y]);
        const double Inertia_J  = Inertia_Ix + Inertia_Iy;

        const double my_gamma    = element->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = neighbor->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        const double k_rot_x = equiv_young * Inertia_Ix / distance;
        const double k_rot_y = equiv_young * Inertia_Iy / distance;
        const double k_tor   = equiv_shear * Inertia_J  / distance;

        double norm_distance = (element->GetRadius() + neighbor->GetRadius()) / distance; // If spheres are not tangent the Damping coefficient, DeltaRotatedAngle and DeltaAngularVelocity have to be normalized

        const double visc_param_rot_x = 2.0 * equiv_gamma * sqrt(equiv_mass * k_rot_x) * norm_distance;
        const double visc_param_rot_y = 2.0 * equiv_gamma * sqrt(equiv_mass * k_rot_y) * norm_distance;
        const double visc_param_tor   = 2.0 * equiv_gamma * sqrt(equiv_mass * k_tor  ) * norm_distance;

        ElasticLocalRotationalMoment[0] = -k_rot_x * LocalDeltaRotatedAngle[0] * norm_distance;
        ElasticLocalRotationalMoment[1] = -k_rot_y * LocalDeltaRotatedAngle[1] * norm_distance;
        ElasticLocalRotationalMoment[2] = -k_tor   * LocalDeltaRotatedAngle[2] * norm_distance;

        ViscoLocalRotationalMoment[0] = -visc_param_rot_x * LocalDeltaAngularVelocity[0] * norm_distance;
        ViscoLocalRotationalMoment[1] = -visc_param_rot_y * LocalDeltaAngularVelocity[1] * norm_distance;
        ViscoLocalRotationalMoment[2] = -visc_param_tor   * LocalDeltaAngularVelocity[2] * norm_distance;

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
