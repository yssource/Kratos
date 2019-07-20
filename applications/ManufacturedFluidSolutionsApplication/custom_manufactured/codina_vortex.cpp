//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "codina_vortex.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

CodinaVortex::CodinaVortex(
    Properties::Pointer pProperties,
    Parameters::Pointer pParameters)
     : ManufacturedSolution(pProperties, pParameters)
{
    // Getting the fluid properties
    mDensity = (*mpProperties)[DENSITY];
    mInvDensity = 1 / mDensity;
    mKinematicViscosity = mInvDensity * (*mpProperties)[DYNAMIC_VISCOSITY];

    // Getting the manufactured settings
    Parameters default_parameters( R"(
    {
        "velocity"    : 1.0,
        "frequency"   : 1.0,
        "damping"     : 1.0
    })");

    mpParameters->ValidateAndAssignDefaults(default_parameters);

    double aux = 16 * std::pow(3, 3/2);
    mLength = 1.0;
    mVelocity = (*mpParameters)["velocity"].GetDouble();
    mOmega = 2 * M_PI * (*mpParameters)["frequency"].GetDouble();
    mDamp = (*mpParameters)["damping"].GetDouble();
    mA = std::sqrt(aux * mVelocity);
}


bool CodinaVortex::IsInsideDomain(array_1d<double, 3>& rCoords)
{
    bool is_inside = true;
    if (rCoords[0] < 0.0 || rCoords[0] > mLength) {is_inside = false;}
    if (rCoords[1] < 0.0 || rCoords[1] > mLength) {is_inside = false;}
    return is_inside;
}


// Velocity components definition
double CodinaVortex::U1(array_1d<double, 3>& rCoords, double& rTime)
{
    return F(rCoords[0]) * DF(rCoords[1]) * G(rTime);
}

double CodinaVortex::U2(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DF(rCoords[0]) * F(rCoords[1]) * G(rTime);
}


// Velocity time derivative definition
double CodinaVortex::DU1_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    return F(rCoords[0]) * DF(rCoords[1]) * DG(rTime);
}

double CodinaVortex::DU2_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DF(rCoords[0]) * F(rCoords[1]) * DG(rTime);
}


// Velocity first derivatives definition
double CodinaVortex::DU1_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return DF(rCoords[0]) * DF(rCoords[1]) * G(rTime);
}

double CodinaVortex::DU1_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return F(rCoords[0]) * DDF(rCoords[1]) * G(rTime);
}

double CodinaVortex::DU2_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DDF(rCoords[0]) * F(rCoords[1]) * G(rTime);
}

double CodinaVortex::DU2_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DF(rCoords[0]) * DF(rCoords[1]) * G(rTime);
}


// Velocity second derivatives definition
double CodinaVortex::DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return DDF(rCoords[0]) * DF(rCoords[1]) * G(rTime);
}

double CodinaVortex::DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return F(rCoords[0]) * DDDF(rCoords[1]) * G(rTime);
}

double CodinaVortex::DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DDDF(rCoords[0]) * F(rCoords[1]) * G(rTime);
}

double CodinaVortex::DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return -DF(rCoords[0]) * DDF(rCoords[1]) * G(rTime);
}


// Auxiliary functions
double CodinaVortex::G(double& rT)
{
    return std::cos(mOmega * rT) * std::exp(-mDamp * rT);
}

double CodinaVortex::DG(double& rT)
{
    return -mOmega * std::sin(mOmega * rT) * std::exp(-mDamp * rT) - mDamp * std::cos(mOmega * rT) * std::exp(-mDamp * rT);
}

double CodinaVortex::F(double& rX)
{
    return mA * std::pow(rX, 2) * std::pow(mLength - rX, 2);
}

double CodinaVortex::DF(double& rX)
{
    return 2.0 * mA * rX * std::pow(mLength - rX, 2) - 2.0 * mA * std::pow(rX, 2) * (mLength - rX);
}

double CodinaVortex::DDF(double& rX)
{
    return 2.0 * mA * (std::pow(mLength, 2) - 6.0 * mLength * rX + 6 * std::pow(rX, 2));
}

double CodinaVortex::DDDF(double& rX)
{
    return 12.0 * mA * (2.0 * rX - mLength);
}


}  // namespace Kratos.


