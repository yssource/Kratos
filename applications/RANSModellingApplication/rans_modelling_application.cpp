//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "rans_modelling_application.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
KratosRANSModellingApplication::KratosRANSModellingApplication()
    : KratosApplication("RANSModellingApplication"),
      mRANSEVMK2D(0,
                  Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                      Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMK3D(0,
                  Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                      Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMEPSILON2D(0,
                        Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                            Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMEPSILON3D(0,
                        Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                            Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMKAdjoint2D(0,
                         Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                             Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMKAdjoint3D(0,
                         Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                             Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMEpsilonAdjoint2D(0,
                               Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                   Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMEpsilonAdjoint3D(0,
                               Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                   Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMKEpsilonVMSAdjoint2D(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMKEpsilonVMSAdjoint3D(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mRANSEVMMonolithicKEpsilonVMSAdjoint2D(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRANSEVMMonolithicKEpsilonVMSAdjoint3D(
          0,
          Element::GeometryType::Pointer(
              new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4))))
{
}

void KratosRANSModellingApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosRANSModellingApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_2)
    KRATOS_REGISTER_VARIABLE(IS_CO_SOLVING_PROCESS_ACTIVE)
    KRATOS_REGISTER_VARIABLE(OLD_CONVERGENCE_VARIABLE)
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_1)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_2)
    KRATOS_REGISTER_VARIABLE(WALL_SMOOTHNESS_BETA)
    KRATOS_REGISTER_VARIABLE(WALL_VON_KARMAN)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C_MU)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MIN)
    KRATOS_REGISTER_VARIABLE(TURBULENT_VISCOSITY_MAX)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA)

    // Register adjoint variables
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS_VELOCITY_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_PRESSURE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_ACCELERATION_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_KINETIC_ENERGY_RATE_PARTIAL_DERIVATIVE)
    KRATOS_REGISTER_VARIABLE(RANS_TURBULENT_ENERGY_DISSIPATION_RATE_2_PARTIAL_DERIVATIVE)

    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_1)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_2)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_1_ADJOINT_3)
    KRATOS_REGISTER_VARIABLE(RANS_AUX_ADJOINT_SCALAR_1 )

    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_1)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_2)
    KRATOS_REGISTER_VARIABLE(RANS_SCALAR_2_ADJOINT_3)
    KRATOS_REGISTER_VARIABLE(RANS_AUX_ADJOINT_SCALAR_2 )

    // Register Elements
    KRATOS_REGISTER_ELEMENT("RANSEVMK2D3N", mRANSEVMK2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMK3D4N", mRANSEVMK3D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON2D3N", mRANSEVMEPSILON2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON3D4N", mRANSEVMEPSILON3D);

    KRATOS_REGISTER_ELEMENT("RANSEVMKAdjoint2D3N", mRANSEVMKAdjoint2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMKAdjoint3D4N", mRANSEVMKAdjoint3D);

    KRATOS_REGISTER_ELEMENT("RANSEVMEpsilonAdjoint2D3N", mRANSEVMEpsilonAdjoint2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEpsilonAdjoint3D4N", mRANSEVMEpsilonAdjoint3D);

    KRATOS_REGISTER_ELEMENT("RANSEVMKEpsilonVMSAdjoint2D3N", mRANSEVMKEpsilonVMSAdjoint2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMKEpsilonVMSAdjoint3D4N", mRANSEVMKEpsilonVMSAdjoint3D);

    KRATOS_REGISTER_ELEMENT("RANSEVMMonolithicKEpsilonVMSAdjoint2D",
                            mRANSEVMMonolithicKEpsilonVMSAdjoint2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMMonolithicKEpsilonVMSAdjoint3D",
                            mRANSEVMMonolithicKEpsilonVMSAdjoint3D);
}
} // namespace Kratos.
