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
    mInvDensity = 1 / mDensity;
    mKinematicViscosity = mInvDensity * (*mpProperties)[DYNAMIC_VISCOSITY];

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


double EcaFlow::Reynolds()
{
    return mVelocity * mLength / mKinematicViscosity;
}


double EcaFlow::Strouhal()
{
    double frequency = 0.5 * mOmega / M_PI;
    return frequency * mLength / mVelocity;
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
    double x = rCoords[0];
    double y = rCoords[1];
    return mA*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) + mVelocity*std::erf(mSigma*y/x);
}

double EcaFlow::U2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 2*M_PI*mA*mVelocity*(-(mB*y + 1)*std::exp(-mB*y) + 1)*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::cos(M_PI*(-2*x + 1))/std::pow(mB, 2) + mVelocity*(1 - std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2)))/(std::sqrt(M_PI)*mSigma);
}


// Velocity time derivative definition
double EcaFlow::DU1_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return mA*mVelocity*y*((mIsPeriodic) ? (
   mOmega*std::sin(mOmega*rTime)
)
: (
   2.5*std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1));
}

double EcaFlow::DU2_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 2*M_PI*mA*mVelocity*(-(mB*y + 1)*std::exp(-mB*y) + 1)*((mIsPeriodic) ? (
   mOmega*std::sin(mOmega*rTime)
)
: (
   2.5*std::exp(-2.5*rTime)
))*std::cos(M_PI*(-2*x + 1))/std::pow(mB, 2);
}


// Velocity first derivatives definition
double EcaFlow::DU1_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return -2*M_PI*mA*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::cos(M_PI*(-2*x + 1)) - 2*mSigma*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 2));
}

double EcaFlow::DU1_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return -mA*mB*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) + mA*mVelocity*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) + 2*mSigma*mVelocity*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*x);
}

double EcaFlow::DU2_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 4*std::pow(M_PI, 2)*mA*mVelocity*(-(mB*y + 1)*std::exp(-mB*y) + 1)*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::sin(M_PI*(-2*x + 1))/std::pow(mB, 2) - 2*mSigma*mVelocity*std::pow(y, 2)*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 3));
}

double EcaFlow::DU2_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 2*M_PI*mA*mVelocity*(-mB*(-mB*y - 1)*std::exp(-mB*y) - mB*std::exp(-mB*y))*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::cos(M_PI*(-2*x + 1))/std::pow(mB, 2) + 2*mSigma*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 2));
}


// Velocity second derivatives definition
double EcaFlow::DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return -4*std::pow(M_PI, 2)*mA*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 4*std::pow(mSigma, 3)*mVelocity*std::pow(y, 3)*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 5)) + 4*mSigma*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 3));
}

double EcaFlow::DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return mA*std::pow(mB, 2)*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 2*mA*mB*mVelocity*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 4*std::pow(mSigma, 3)*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 3));
}

double EcaFlow::DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return -4*std::pow(M_PI, 2)*mA*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 4*std::pow(mSigma, 3)*mVelocity*std::pow(y, 3)*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 5)) + 4*mSigma*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 3));
}

double EcaFlow::DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return mA*std::pow(mB, 2)*mVelocity*y*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 2*mA*mB*mVelocity*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::exp(-mB*y)*std::sin(M_PI*(-2*x + 1)) - 4*std::pow(mSigma, 3)*mVelocity*y*std::exp(-std::pow(mSigma, 2)*std::pow(y, 2)/std::pow(x, 2))/(std::sqrt(M_PI)*std::pow(x, 3));
}


// Pressure and derivatives definition
double EcaFlow::P(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 0.050000000000000003*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::sin(M_PI*(8.0*std::pow(y, 3) - 6.0*std::pow(y, 2) + 0.5))*std::sin(M_PI*(0.5*std::pow(4.0*x - 3.0, 4) - 1.0*std::pow(4.0*x - 3.0, 2) + 0.5)) + 50.0*std::log(0.5*std::pow(y, 2) + 0.875)*std::log(0.33333333333333331*std::pow(x, 3) - 0.75*std::pow(x, 2) + 0.5*x + 0.91666666666666663);
}

double EcaFlow::DP_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 0.050000000000000003*M_PI*(-32.0*x + 8.0*std::pow(4.0*x - 3.0, 3) + 24.0)*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::sin(M_PI*(8.0*std::pow(y, 3) - 6.0*std::pow(y, 2) + 0.5))*std::cos(M_PI*(0.5*std::pow(4.0*x - 3.0, 4) - 1.0*std::pow(4.0*x - 3.0, 2) + 0.5)) + 50.0*(1.0*std::pow(x, 2) - 1.5*x + 0.5)*std::log(0.5*std::pow(y, 2) + 0.875)/(0.33333333333333331*std::pow(x, 3) - 0.75*std::pow(x, 2) + 0.5*x + 0.91666666666666663);
}

double EcaFlow::DP_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return 50.0*y*std::log(0.33333333333333331*std::pow(x, 3) - 0.75*std::pow(x, 2) + 0.5*x + 0.91666666666666663)/(0.5*std::pow(y, 2) + 0.875) + 0.050000000000000003*M_PI*(24.0*std::pow(y, 2) - 12.0*y)*((mIsPeriodic) ? (
   -std::cos(mOmega*rTime) + 1.0
)
: (
   1.0 - std::exp(-2.5*rTime)
))*std::sin(M_PI*(0.5*std::pow(4.0*x - 3.0, 4) - 1.0*std::pow(4.0*x - 3.0, 2) + 0.5))*std::cos(M_PI*(8.0*std::pow(y, 3) - 6.0*std::pow(y, 2) + 0.5));
}


}  // namespace Kratos.


