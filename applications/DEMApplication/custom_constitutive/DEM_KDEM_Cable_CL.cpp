// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Cable_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Cable::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Cable(*this));
        return p_clone;
    }

    void DEM_KDEM_Cable::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Cable to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_Cable::Check(Properties::Pointer pProp) const {
        DEMContinuumConstitutiveLaw::Check(pProp);

        if(!pProp->Has(ROTATIONAL_MOMENT_COEFFICIENT)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable ROTATIONAL_MOMENT_COEFFICIENT should be present in the properties when using DEM_KDEM. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(ROTATIONAL_MOMENT_COEFFICIENT) = 1.0;
        }
    }

    void DEM_KDEM_Cable::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = 0.5 * radius_sum;
        calculation_area = Globals::Pi * equiv_radius * equiv_radius;
        KRATOS_CATCH("")
    }

    void DEM_KDEM_Cable::GetContactArea(const double radius, const double other_radius, const Vector& vector_of_initial_areas, const int neighbour_position, double& calculation_area) {
        if (vector_of_initial_areas.size()) calculation_area = vector_of_initial_areas[neighbour_position];
        else CalculateContactArea(radius, other_radius, calculation_area);
    }

    void DEM_KDEM_Cable::CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                                    double& equiv_visco_damp_coeff_tangential,
                                                    SphericContinuumParticle* element1,
                                                    SphericContinuumParticle* element2,
                                                    const double kn_el,
                                                    const double kt_el) {

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
        const double equiv_mass = 0.5 * (my_mass + other_mass);
        // const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect) = element1->GetGeometry()[0].Coordinates() - element2->GetGeometry()[0].Coordinates();
        const double distance = DEM_MODULUS_3(other_to_me_vect);
        const double norm_distance = (element1->GetRadius() + element2->GetRadius()) / distance; // If spheres are not tangent the Damping coefficient has to be normalized

        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn_el) * norm_distance;
        equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * kt_el) * norm_distance;
    }

    void DEM_KDEM_Cable::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                          SphericContinuumParticle* neighbor,
                                                          double equiv_young,
                                                          double distance,
                                                          double calculation_area,
                                                          double LocalCoordSystem[3][3],
                                                          double ElasticLocalRotationalMoment[3],
                                                          double ViscoLocalRotationalMoment[3],
                                                          double equiv_poisson,
                                                          double indentation) {

        double LocalDeltaRotatedAngle[3]    = {0.0};
        double LocalDeltaAngularVelocity[3] = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3> GlobalDeltaAngularVelocity;
        noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);
        //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mContactMoment, LocalRotationalMoment);

        const double equivalent_radius = sqrt(calculation_area / Globals::Pi);

        const double element_mass  = element->GetMass();
        const double neighbor_mass = neighbor->GetMass();
        // const double equiv_mass    = element_mass * neighbor_mass / (element_mass + neighbor_mass);
        const double equiv_mass    = 0.5 * (element_mass + neighbor_mass);

        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double Inertia_I     = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
        const double Inertia_J     = 2.0 * Inertia_I; // This is the polar inertia

        const double my_gamma    = element->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = neighbor->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        const double k_rot = equiv_young * Inertia_I / distance;
        const double k_tor = equiv_shear * Inertia_J / distance;

        const double Moment_of_inertia_I = 0.25 * equiv_mass * pow(equivalent_radius, 2.0) + 0.333333333333333 * equiv_mass * pow(distance, 2.0);
        const double Moment_of_inertia_J = 0.5 * equiv_mass * equivalent_radius * equivalent_radius;

        const double norm_distance = (element->GetRadius() + neighbor->GetRadius()) / distance; // If spheres are not tangent the Damping coefficient, DeltaRotatedAngle and DeltaAngularVelocity have to be normalized

        //Viscous parameter taken from Olmedo et al., 'Discrete element model of the dynamic response of fresh wood stems to impact'
        array_1d<double, 3> visc_param;
        const double visc_param_rot = 2.0 * equiv_gamma * sqrt(Moment_of_inertia_I * k_rot) * norm_distance;
        const double visc_param_tor = 2.0 * equiv_gamma * sqrt(Moment_of_inertia_J * k_tor) * norm_distance;
        // const double visc_param_rot = 2.0 * equiv_gamma * sqrt(Inertia_I * k_rot) * norm_distance;
        // const double visc_param_tor = 2.0 * equiv_gamma * sqrt(Inertia_J * k_tor) * norm_distance;

        ElasticLocalRotationalMoment[0] = -k_rot * LocalDeltaRotatedAngle[0] * norm_distance;
        ElasticLocalRotationalMoment[1] = -k_rot * LocalDeltaRotatedAngle[1] * norm_distance;
        ElasticLocalRotationalMoment[2] = -k_tor * LocalDeltaRotatedAngle[2] * norm_distance;

        ViscoLocalRotationalMoment[0] = -visc_param_rot * LocalDeltaAngularVelocity[0] * norm_distance;
        ViscoLocalRotationalMoment[1] = -visc_param_rot * LocalDeltaAngularVelocity[1] * norm_distance;
        ViscoLocalRotationalMoment[2] = -visc_param_tor * LocalDeltaAngularVelocity[2] * norm_distance;

        // TODO: Judge if the rotation spring is broken or not
        /*
        double ForceN  = LocalElasticContactForce[2];
        double ForceS  = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
        double MomentN = LocalRotaSpringMoment[2];
        // bending stress and axial stress add together, use edge of the bar will failure first
        double TensiMax = -ForceN / calculation_area + MomentS / Inertia_I * equiv_radius;
        double ShearMax =  ForceS / calculation_area + fabs(MomentN) / Inertia_J * equiv_radius;
        if (TensiMax > equiv_tension || ShearMax > equiv_cohesion) {
            mRotaSpringFailureType[i_neighbor_count] = 1;
            LocalRotaSpringMoment[0] = LocalRotaSpringMoment[1] = LocalRotaSpringMoment[2] = 0.0;
            //LocalRotaSpringMoment[1] = 0.0;
            //LocalRotaSpringMoment[2] = 0.0;
        }
        */
        //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotationalMoment, mContactMoment);
    }//ComputeParticleRotationalMoments

} // namespace Kratos
