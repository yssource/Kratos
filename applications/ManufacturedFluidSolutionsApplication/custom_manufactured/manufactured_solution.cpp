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
#include "manufactured_solution.h"


namespace Kratos
{

array_1d<double, 3> ManufacturedSolution::BodyForce(array_1d<double, 3>& rCoords, double& rTime)
{
    KRATOS_ERROR_IF_NOT(IsInsideDomain(rCoords)) << "The given coordinates are outside the domain. Coordinates are : " << rCoords << " Please, check the definition of the manufactured solution" << std::endl;
    auto time_der = TimeDerivative(rCoords, rTime);
    auto convective = ConvectiveTerm(rCoords, rTime);
    auto viscous = ViscousTerm(rCoords, rTime);
    auto press_grad = mInvDensity * PressureGradient(rCoords, rTime);
    return time_der + convective + press_grad - viscous;
}


array_1d<double, 3> ManufacturedSolution::Velocity(array_1d<double, 3>& rCoords, double& rTime)
{
    KRATOS_ERROR_IF_NOT(IsInsideDomain(rCoords)) << "The given coordinates are outside the domain. Coordinates are : " << rCoords << " Please, check the definition of the manufactured solution" << std::endl;
    array_1d<double, 3> velocity;
    velocity[0] = U1(rCoords, rTime);
    velocity[1] = U2(rCoords, rTime);
    velocity[2] = U3(rCoords, rTime);
    return velocity;
}


double ManufacturedSolution::Pressure(array_1d<double, 3>& rCoords, double& rTime)
{
    KRATOS_ERROR_IF_NOT(IsInsideDomain(rCoords)) << "The given coordinates are outside the domain. Coordinates are : " << rCoords << " Please, check the definition of the manufactured solution" << std::endl;
    return P(rCoords, rTime);
}


array_1d<double, 3> ManufacturedSolution::ConvectiveTerm(array_1d<double, 3>& rCoords, double& rTime)
{
    auto vel = Velocity(rCoords, rTime);
    auto grad = VelocityGradient(rCoords, rTime);
    return prod(vel, grad);
}


array_1d<double, 3> ManufacturedSolution::ViscousTerm(array_1d<double, 3>& rCoords, double& rTime)
{
    return mKinematicViscosity * VelocityLaplacian(rCoords, rTime);
}


array_1d<double, 3> ManufacturedSolution::TimeDerivative(array_1d<double, 3>& rCoords, double& rTime)
{
    array_1d<double, 3> time_der;
    time_der[0] = DU1_DT(rCoords, rTime);
    time_der[1] = DU2_DT(rCoords, rTime);
    time_der[2] = DU3_DT(rCoords, rTime);
    return time_der;
}


BoundedMatrix<double, 3, 3> ManufacturedSolution::VelocityGradient(array_1d<double, 3>& rCoords, double& rTime)
{
    BoundedMatrix<double, 3, 3> grad;
    grad(0,0) = DU1_DX1(rCoords, rTime);
    grad(1,0) = DU1_DX2(rCoords, rTime);
    grad(2,0) = DU1_DX3(rCoords, rTime);
    grad(0,1) = DU2_DX1(rCoords, rTime);
    grad(1,1) = DU2_DX2(rCoords, rTime);
    grad(2,1) = DU2_DX3(rCoords, rTime);
    grad(0,2) = DU3_DX1(rCoords, rTime);
    grad(1,2) = DU3_DX2(rCoords, rTime);
    grad(2,2) = DU3_DX3(rCoords, rTime);
    return grad;
}


array_1d<double, 3> ManufacturedSolution::VelocityLaplacian(array_1d<double, 3>& rCoords, double& rTime)
{
    array_1d<double, 3> laplacian;
    laplacian[0] = DDU1_DX11(rCoords, rTime) + DDU1_DX22(rCoords, rTime) + DDU1_DX33(rCoords, rTime);
    laplacian[1] = DDU2_DX11(rCoords, rTime) + DDU2_DX22(rCoords, rTime) + DDU2_DX33(rCoords, rTime);
    laplacian[2] = DDU3_DX11(rCoords, rTime) + DDU3_DX22(rCoords, rTime) + DDU3_DX33(rCoords, rTime);
    return laplacian;
}


array_1d<double, 3> ManufacturedSolution::PressureGradient(array_1d<double, 3>& rCoords, double& rTime)
{
    array_1d<double, 3> grad;
    grad[0] = DP_DX1(rCoords, rTime);
    grad[1] = DP_DX2(rCoords, rTime);
    grad[2] = DP_DX3(rCoords, rTime);
    return grad;
}


Parameters& ManufacturedSolution::GetParameters()
{
    return *mpParameters;
}


Properties& ManufacturedSolution::GetProperties()
{
    return *mpProperties;
}

}  // namespace Kratos.
