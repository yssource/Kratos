//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "rans_constitutive_laws_application.h"
#include "rans_constitutive_laws_application_variables.h"


namespace Kratos {

KratosRANSConstitutiveLawsApplication::KratosRANSConstitutiveLawsApplication():
    KratosApplication("RANSConstitutiveLawsApplication"),
    mRANSEVMK2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mRANSEVMK3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mRANSEVMEPSILON2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mRANSEVMEPSILON3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mEpsilonWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mEpsilonWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
    {}

void KratosRANSConstitutiveLawsApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosRANSConstitutiveLawsApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_VARIABLE( RANS_MODELLING_PROCESS_STEP )
    KRATOS_REGISTER_VARIABLE( RANS_Y_PLUS )

    // Turbulence model constants
    KRATOS_REGISTER_VARIABLE( WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_VARIABLE( WALL_VON_KARMAN )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C_MU )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C1 )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_RANS_C2 )
    KRATOS_REGISTER_VARIABLE( TURBULENT_MIXING_LENGTH )
    KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY_FRACTION )
    KRATOS_REGISTER_VARIABLE( TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_REGISTER_VARIABLE( TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_VARIABLE( PARENT_ELEMENT )
    KRATOS_REGISTER_VARIABLE( RANS_VELOCITY_STEP )
    KRATOS_REGISTER_VARIABLE( GAUSS_POINT_INDICES)

    // Register Elements
    KRATOS_REGISTER_ELEMENT("RANSEVMK2D3N",mRANSEVMK2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMK3D4N",mRANSEVMK3D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON2D3N",mRANSEVMEPSILON2D);
    KRATOS_REGISTER_ELEMENT("RANSEVMEPSILON3D4N",mRANSEVMEPSILON3D);

    KRATOS_REGISTER_CONDITION("EpsilonWallCondition2D", mEpsilonWallCondition2D);
    KRATOS_REGISTER_CONDITION("EpsilonWallCondition3D", mEpsilonWallCondition2D);
}
}  // namespace Kratos.
