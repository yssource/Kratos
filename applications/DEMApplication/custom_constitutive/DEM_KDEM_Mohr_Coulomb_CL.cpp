// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Mohr_Coulomb_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Mohr_Coulomb::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Mohr_Coulomb(*this));
        return p_clone;
    }

    void DEM_KDEM_Mohr_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Mohr_Coulomb to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_Mohr_Coulomb::Check(Properties::Pointer pProp) const {
        DEM_KDEM::Check(pProp);

        if(!pProp->Has(INTERNAL_COHESION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable INTERNAL_COHESION should be present in the properties when using DEM_KDEM_Mohr_Coulomb. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(INTERNAL_COHESION) = 0.0;
        }
        if(!pProp->Has(INTERNAL_FRICTION_ANGLE)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable INTERNAL_FRICTION_ANGLE should be present in the properties when using DEM_KDEM_Mohr_Coulomb. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(INTERNAL_FRICTION_ANGLE) = 0.0;
        }
    }

    double DEM_KDEM_Mohr_Coulomb::LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);

        GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);

        const double radius_sum = my_radius + other_radius;
        const double initial_delta = element1->GetInitialDelta(i);
        const double initial_dist = radius_sum - initial_delta;
        const double kn_el = equiv_young * calculation_area / initial_dist;

        const double max_normal_force = mohr_coulomb_c * calculation_area;
        const double u1 = max_normal_force / kn_el;
        return u1;
    }

    void DEM_KDEM_Mohr_Coulomb::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps) {

        KRATOS_TRY

        if (indentation >= 0.0) { //COMPRESSION
            LocalElasticContactForce[2] = kn_el * indentation;
        }
        else { //tension
            int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
            if (failure_type == 0) {
                LocalElasticContactForce[2] = kn_el * indentation;
            }
            else {
                LocalElasticContactForce[2] = 0.0;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Mohr_Coulomb::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
            double& contact_tau,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            DenseVector<int>& search_control_vector,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential

        double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                  + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        if (failure_type == 0) { // This means it has not broken
            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) { //TODO: use this only for intact bonds (not broken))
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }
            contact_tau = ShearForceNow / calculation_area;
            contact_sigma = LocalElasticContactForce[2] / calculation_area;
        }
        else {
            const double equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Mohr_Coulomb::CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2){

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (failure_type == 0) {
            int use_element_tensor_1, use_element_tensor_2 = 0;

            if(element1->IsSkin() && element2->IsSkin()) {
                if (element1->IsNot(DEMFlags::COPIED_STRESS_TENSOR) &&  element1->IsNot(DEMFlags::COPIED_STRESS_TENSOR)) {
                    if(element2->IsNot(DEMFlags::COPIED_STRESS_TENSOR) &&  element2->IsNot(DEMFlags::COPIED_STRESS_TENSOR)) {
                        use_element_tensor_1 = 0;
                        use_element_tensor_2 = 0;
                    }
                    else {
                        use_element_tensor_1 = 0;
                        use_element_tensor_2 = 1;
                    }
                }
                else {
                    if (element2->IsNot(DEMFlags::COPIED_STRESS_TENSOR) &&  element2->IsNot(DEMFlags::COPIED_STRESS_TENSOR)) {
                        use_element_tensor_1 = 1;
                        use_element_tensor_2 = 0;
                    }
                    else {
                        use_element_tensor_1 = 1;
                        use_element_tensor_2 = 1;
                    }
                }
            }
            else {
                if (element1->IsSkin()) {
                    use_element_tensor_1 = 0;
                    use_element_tensor_2 = 1;
                }
                else {
                    use_element_tensor_1 = 1;
                    use_element_tensor_2 = 0;
                }
            }

            CheckMohrCoulombFailure(i_neighbour_count, element1, element2, use_element_tensor_1, use_element_tensor_2);
        }
    }

    void DEM_KDEM_Mohr_Coulomb::CheckMohrCoulombFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2, int use_element_tensor_1, int use_element_tensor_2) {
        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();

        if ( !use_element_tensor_1 && !use_element_tensor_2) {
            const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);

        }
        else {
            Matrix average_stress_tensor = ZeroMatrix(3,3);
            double factor = 1.0 / (use_element_tensor_1 + use_element_tensor_2);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    average_stress_tensor(i,j) = factor * (use_element_tensor_1 *(*(element1->mSymmStressTensor))(i,j) + use_element_tensor_2 *(*(element2->mSymmStressTensor))(i,j));
                }
            }

            Vector principal_stresses(3);
            noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);

            Properties& element1_props = element1->GetProperties();
            Properties& element2_props = element2->GetProperties();

            const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);
            const double mohr_coulomb_phi = 0.5 * (element1_props[INTERNAL_FRICTION_ANGLE] + element2_props[INTERNAL_FRICTION_ANGLE]);
            const double mohr_coulomb_phi_in_radians = mohr_coulomb_phi * Globals::Pi / 180.0;
            const double sinphi = std::sin(mohr_coulomb_phi_in_radians);
            const double cosphi = std::cos(mohr_coulomb_phi_in_radians);

            const double max_stress = *std::max_element(principal_stresses.begin(), principal_stresses.end());
            const double min_stress = *std::min_element(principal_stresses.begin(), principal_stresses.end());
            const double function_value = (max_stress - min_stress) + (max_stress + min_stress) * sinphi - 2.0 * mohr_coulomb_c * cosphi;

            if(function_value > 0) {
                int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
                failure_type = 4;
            }
        }
    }

    bool DEM_KDEM_Mohr_Coulomb::CheckRequirementsOfStressTensor() {

        return true;

    }

} // namespace Kratos
