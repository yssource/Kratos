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
        "is_periodic" : true,
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
    return -2.0 / std::sqrt(M_PI) * mSigma * rCoords[1] / std::pow(rCoords[0], 2) * ExpEta2(rCoords)
        - 2.0 * M_PI * mA * rCoords[1] * ExpBY(rCoords) * CosFunc(rCoords) * F(rTime);
}

double EcaFlow::DU1_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2.0 / std::sqrt(M_PI) * mSigma / rCoords[0] * ExpEta2(rCoords)
        + mA * ExpBY(rCoords) * SinFunc(rCoords) * (1.0 - mB * rCoords[1]) * F(rTime);
}

double EcaFlow::DU2_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return -2.0 / std::sqrt(M_PI) * mSigma * std::pow(rCoords[1], 2) / std::pow(rCoords[0], 3) * ExpEta2(rCoords)
        + 4.0 * mA * std::pow(M_PI, 2) / std::pow(mB, 2) * SinFunc(rCoords) * (1.0 - ExpBY(rCoords) * (mB * rCoords[1] + 1.0)) * F(rTime);
}

double EcaFlow::DU2_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2.0 / std::sqrt(M_PI) * mSigma * rCoords[1] / std::pow(rCoords[0], 2) * ExpEta2(rCoords)
        + 2.0 * mA * M_PI * rCoords[1] * ExpBY(rCoords) * CosFunc(rCoords) * F(rTime);
}


// Velocity second derivatives definition
double EcaFlow::DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return 4.0 / std::sqrt(M_PI) * Eta(rCoords) / std::pow(rCoords[0], 2) * ExpEta2(rCoords) * (1.0 - std::pow(Eta(rCoords), 2))
        -4.0 * mA * std::pow(M_PI, 2) * rCoords[1] * ExpBY(rCoords) * SinFunc(rCoords) * F(rTime);
}

double EcaFlow::DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return -4.0 / std::sqrt(M_PI) * std::pow(mSigma / rCoords[0], 2) * Eta(rCoords) * ExpEta2(rCoords)
        + mA * mB * ExpBY(rCoords) * SinFunc(rCoords) * (rCoords[1] * mB - 2.0) * F(rTime);
}

double EcaFlow::DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2.0 / std::sqrt(M_PI) * mSigma * std::pow(rCoords[1], 2) / std::pow(rCoords[0], 3) * ExpEta2(rCoords) * (3.0 - 2.0 * std::pow(Eta(rCoords), 2))
        + 8.0 * mA * std::pow(M_PI, 3) / std::pow(mB, 2) * CosFunc(rCoords) * (ExpBY(rCoords) * (mB * rCoords[1] + 1.0) - 1.0) * F(rTime);
}

double EcaFlow::DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    return 2.0 / std::sqrt(M_PI) * mSigma / std::pow(rCoords[0], 2) * ExpEta2(rCoords) * (1.0 - 2.0 * std::pow(Eta(rCoords), 2))
        + 2.0 * mA * M_PI * ExpBY(rCoords) * CosFunc(rCoords) * (1 - mB * rCoords[1]) * F(rTime);
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
    return mSigma * rCoords[1] / rCoords[0];
}

double EcaFlow::SinFunc(array_1d<double, 3>& rCoords)
{
    return std::sin((1.0 - 2.0 * rCoords[0]) * M_PI);
}

double EcaFlow::CosFunc(array_1d<double, 3>& rCoords)
{
    return std::cos((1.0 - 2.0 * rCoords[0]) * M_PI);
}

double EcaFlow::ExpBY(array_1d<double, 3>& rCoords)
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
    return mA * rCoords[1] * ExpBY(rCoords) * SinFunc(rCoords);
}

double EcaFlow::UB2(array_1d<double, 3>& rCoords)
{
    return 2 * mA * M_PI / std::pow(mB, 2) * CosFunc(rCoords) * (1.0 - ExpBY(rCoords) * (mB * rCoords[1] + 1.0));
}


// Pressure and derivatives definition
double EcaFlow::P(array_1d<double, 3>& rCoords, double& rTime)
{
    return 50.0 * std::log(PA1(rCoords)) * std::log(PA2(rCoords)) 
        -0.05 * std::sin(PB1(rCoords) * M_PI / 2.0) * std::sin(PB2(rCoords) * M_PI / 2.0) * F(rTime);
}

double EcaFlow::DP_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    return 50.0 * DPA1(rCoords) * std::log(PA2(rCoords)) / PA1(rCoords)
        -0.05 * DPB1(rCoords) * M_PI / 2.0 * std::cos(PB1(rCoords) * M_PI / 2.0) * std::sin(PB2(rCoords) * M_PI / 2.0) * F(rTime);
}

double EcaFlow::DP_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    return 50.0 * DPA2(rCoords) * std::log(PA1(rCoords)) / PA2(rCoords)
        -0.05 * DPB2(rCoords) * M_PI / 2.0 * std::sin(PB1(rCoords) * M_PI / 2.0) * std::cos(PB2(rCoords) * M_PI / 2.0) * F(rTime);
}


// Auxiliary functions for the pressure
double EcaFlow::PA1(array_1d<double, 3>& rCoords)
{
    double x = rCoords[0];
    return std::pow(x, 3) / 3.0 - 3.0 * std::pow(x, 2) / 4.0 + x / 2.0 + 11.0 / 12.0;
}

double EcaFlow::PA2(array_1d<double, 3>& rCoords)
{
    double y = rCoords[1];
    return std::pow(y, 2) / 2 + 7.0 / 8.0;
}

double EcaFlow::PB1(array_1d<double, 3>& rCoords)
{
    double x = rCoords[0];
    return std::pow(4.0 * x - 3.0, 4) - 2.0 * std::pow(4.0 * x - 3.0, 2) + 1.0;
}

double EcaFlow::PB2(array_1d<double, 3>& rCoords)
{
    double y = rCoords[1];
    return 16.0 * std::pow(y, 3) - 12.0 * std::pow(y, 2) + 1.0;
}


// Derivatives of the auxiliary functions for the pressure
double EcaFlow::DPA1(array_1d<double, 3>& rCoords)
{
    double x = rCoords[0];
    return std::pow(x, 2) - 3.0 * x / 2.0 + 1.0 / 2.0;
}

double EcaFlow::DPA2(array_1d<double, 3>& rCoords)
{
    double y = rCoords[1];
    return y;
}

double EcaFlow::DPB1(array_1d<double, 3>& rCoords)
{
    double x = rCoords[0];
    return 16.0 * std::pow(4.0 * x - 3.0, 3) - 16.0 * std::pow(4.0 * x - 3.0, 2);
}

double EcaFlow::DPB2(array_1d<double, 3>& rCoords)
{
    double y = rCoords[1];
    return 48.0 * std::pow(y, 2) - 24.0 * y;
}

}  // namespace Kratos.


