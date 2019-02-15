//
// Author: Marc Chung To Sang mchungtosang@cimne.upc.edu
//          

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "ion_particle.h"

namespace Kratos {

IonParticle::~IonParticle(){}

IonParticle& IonParticle::operator=(IonParticle const& rOther) {
    
    SphericParticle::operator=(rOther);
        
    mSingleIonCharge = rOther.mSingleIonCharge;
    mDoubleIonCharge = rOther.mDoubleIonCharge;
    mXenonMass = rOther.mXenonMass;

    return *this;
}

void IonParticle::Initialize(const ProcessInfo& r_process_info) {   
    SphericParticle::Initialize(r_process_info);
    double added_mass_coefficient = 1.0;
    SetMass(added_mass_coefficient * GetDensity() * CalculateVolume());
}

void IonParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                             array_1d<double, 3>& additionally_applied_moment,
                             const ProcessInfo& r_current_process_info,
                             const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    array_1d<double, 3> Coulomb_force; 
    CalculateCoulombForce(Coulomb_force);

    array_1d<double, 3> Laplace_force;          //Laplace_force.clear();
    CalculateLaplaceForce(Laplace_force);

    noalias(additionally_applied_force) += Coulomb_force + Laplace_force;

    //Now add the contribution of base class function (gravity or other forces added in upper levels):
    SphericParticle::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);
    KRATOS_CATCH( "" )
}

void IonParticle::CalculateCoulombForce(array_1d<double, 3>& Coulomb_force)
{
    double Density =  this->GetDensity();
    Coulomb_force[0] =  Density * (mSingleIonCharge / mXenonMass) * mExternalElectricField[0];
    Coulomb_force[1] =  Density * (mSingleIonCharge / mXenonMass) * mExternalElectricField[1];
    Coulomb_force[2] =  Density * (mSingleIonCharge / mXenonMass) * mExternalElectricField[2];
}


void IonParticle::CalculateLaplaceForce(array_1d<double, 3>& Laplace_force)
{
    array_1d<double, 3 >& VelocityPreviousStep = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    double Density =  this->GetDensity();
    DEM_SET_TO_CROSS_OF_FIRST_TWO_3(VelocityPreviousStep, mExternalMagneticField, Laplace_force)
    Laplace_force *= Density * (mSingleIonCharge / mXenonMass);
    // Laplace_force[0] = coeff * (VelocityPreviousStep[1] * mExternalMagneticField[2] - VelocityPreviousStep[2] * mExternalMagneticField[1]);
    // Laplace_force[1] = coeff * (VelocityPreviousStep[2] * mExternalMagneticField[0] - VelocityPreviousStep[0] * mExternalMagneticField[2]);
    // Laplace_force[2] = coeff * (VelocityPreviousStep[0] * mExternalMagneticField[1] - VelocityPreviousStep[1] * mExternalMagneticField[0]);
}

void IonParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) 
{
    SphericParticle::MemberDeclarationFirstStep(r_process_info);
}

double IonParticle::GetSingleIonCharge()
{
    return mSingleIonCharge;
}

double IonParticle::GetDoubleIonCharge()
{
    return mDoubleIonCharge;
}

double IonParticle::GetXenonMass()
{
    return mXenonMass;
}

array_1d<double, 3> IonParticle::GetExternalElectricField()
{
    return mExternalElectricField;
}

array_1d<double, 3> IonParticle::GetExternalMagneticField()
{
    return mExternalMagneticField;
}

    
} // namespace Kratos
