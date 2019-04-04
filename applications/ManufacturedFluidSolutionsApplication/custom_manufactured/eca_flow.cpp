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
#include "eca_flow.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

EcaFlow::EcaFlow(
    Properties::Pointer pProperties,
    Parameters::Pointer pParameters)
     : ManufacturedSolution(pProperties, pParameters)
{
    // Getting the fluid properties
    mDensity = (*mpProperties)[DENSITY];
    mDynamicViscosity = (*mpProperties)[DYNAMIC_VISCOSITY];
    mInvDensity = 1 / mDensity;

    // Getting the manufactured settings
    Parameters default_parameters( R"(
    {
        "velocity"    : 1.0,
        "sigma"       : 4.0,
        "a"           : 15,
        "b"           : 20,
        "is_periodic" : True,
        "frequency"   : 1.0
    })");

    mpParameters->ValidateAndAssignDefaults(default_parameters);

    mLength = 1.0;
    mVelocity = (*mpParameters)["velocity"].GetDouble();
    mSigma = (*mpParameters)["sigma"].GetDouble();
    mA = (*mpParameters)["a"].GetDouble();
    mB = (*mpParameters)["b"].GetDouble();
    mIsPeriodic = (*mpParameters)["is_periodic"].GetBool();
    mOmega = 2 * M_PI * (*mpParameters)["frequency"].GetDouble();
}


bool EcaFlow::IsInsideDomain(array_1d<double, 3>& rCoords)
{
    bool is_inside = true;
    if (rCoords[0] < 0.5 || rCoords[0] > mLength) {is_inside = false;}
    if (rCoords[1] < 0.0 || rCoords[1] > 0.5 * mLength) {is_inside = false;}
    return is_inside;
}


// Velocity components definition
double EcaFlow::U1(array_1d<double, 3>& rCoords, double& rTime)
{
    return UA1(rCoords) + UB1(rCoords) * F(rTime);
}

double EcaFlow::U2(array_1d<double, 3>& rCoords, double& rTime)
{
    return UA2(rCoords) + UB2(rCoords) * F(rTime);
}


// Velocity time derivative definition
double EcaFlow::DU1_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    return UB1(rCoords) * DF(rTime);
}

double EcaFlow::DU2_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    return UB2(rCoords) * DF(rTime);
}


// Velocity first derivatives definition
double EcaFlow::DU1_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return -2 / std::sqrt(M_PI) * mSigma * rCoords[1] / std::pow(rCoords[0], 2) * ExpEta2(rCoords)
        - 2 * M_PI * mA * rCoords[1] * ExpBX2(rCoords) * CosFunc(rCoords) * F(rTime);
}

double EcaFlow::DU1_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2 / std::sqrt(M_PI) * mSigma / rCoords[0] * ExpEta2(rCoords)
        + mA * ExpBX2(rCoords) * SinFunc(rCoords) * (1.0 - mB * rCoords[1]) * F(rTime);
}

double EcaFlow::DU2_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return -2 / std::sqrt(M_PI) * mSigma * std::pow(rCoords[1], 2) / std::pow(rCoords[0], 3) * ExpEta2(rCoords)
        + 4 * mA * std::pow(M_PI, 2) / std::pow(mB, 2) * SinFunc(rCoords) * (1.0 - ExpBX2(rCoords) * (mB * rCoords[1] + 1)) * F(rTime);
}

double EcaFlow::DU2_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2 / std::sqrt(M_PI) * mSigma * rCoords[1] / std::pow(rCoords[0], 2) * ExpEta2(rCoords)
        + 2 * mA * M_PI * rCoords[1] * ExpBX2(rCoords) * CosFunc(rCoords) * F(rTime);
}


// Velocity second derivatives definition
double EcaFlow::DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return 0.0;
}

double EcaFlow::DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return 0.0;
}

double EcaFlow::DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return 0.0;
}

double EcaFlow::DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return 0.0;
}


// Auxiliary time dependent functions
double EcaFlow::F(double& rT)
{
    double func;
    if (mIsPeriodic) {func = Fp(rT);} else {func = Fe(rT);}
    return func;
}

double EcaFlow::DF(double& rT)
{
    double der;
    if (mIsPeriodic) {der = DFp(rT);} else {der = DFe(rT);}
    return der;
}

double EcaFlow::Fe(double& rT)
{
    return 1.0 - std::exp(-2.5 * rT);
}

double EcaFlow::DFe(double& rT)
{
    return 2.5 * std::exp(-2.5 * rT);
}

double EcaFlow::Fp(double& rT)
{
    return 1.0 - std::cos(mOmega * rT);
}

double EcaFlow::DFp(double& rT)
{
    return mOmega * std::sin(mOmega * rT);
}


// Auxiliary functions
double EcaFlow::Eta(array_1d<double, 3>& rCoords)
{
    return mSigma * rCoords[0] / rCoords[1];
}

double EcaFlow::SinFunc(array_1d<double, 3>& rCoords)
{
    return std::sin((1 - 2 * rCoords[0]) * M_PI);
}

double EcaFlow::CosFunc(array_1d<double, 3>& rCoords)
{
    return std::cos((1 - 2 * rCoords[0]) * M_PI);
}

double EcaFlow::ExpBX2(array_1d<double, 3>& rCoords)
{
    return std::exp(-mB * rCoords[1]);
}

double EcaFlow::ExpEta2(array_1d<double, 3>& rCoords)
{
    return std::exp(-std::pow(Eta(rCoords), 2));
}


// Auxiliary functions for the velocity definition
double EcaFlow::UA1(array_1d<double, 3>& rCoords)
{
    return std::erf(Eta(rCoords));
}

double EcaFlow::UA2(array_1d<double, 3>& rCoords)
{
    return 1 / (mSigma * std::sqrt(M_PI)) * (1.0 - ExpEta2(rCoords));
}

double EcaFlow::UB1(array_1d<double, 3>& rCoords)
{
    return mA * rCoords[1] * ExpBX2(rCoords) * SinFunc(rCoords);
}

double EcaFlow::UB2(array_1d<double, 3>& rCoords)
{
    return 2 * mA * M_PI / std::pow(mB, 2) * CosFunc(rCoords) * (1.0 - ExpBX2(rCoords) * (mB * rCoords[1] + 1.0));
}


}  // namespace Kratos.


