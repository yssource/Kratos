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
    mKEpsilon2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mKEpsilon3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mKEpsilonWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mKEpsilonWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) )
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
    KRATOS_REGISTER_VARIABLE( ELEMENT_DERIVATIVES_DOFS_EXTENSION )

    // Register Elements
    KRATOS_REGISTER_ELEMENT("KEPSILON2D3N",mKEpsilon2D);
    KRATOS_REGISTER_ELEMENT("KEPSILON3D4N",mKEpsilon3D);

    KRATOS_REGISTER_CONDITION("KEpsilonWallCondition2D2N", mKEpsilonWallCondition2D);
    KRATOS_REGISTER_CONDITION("KEpsilonWallCondition3D3N", mKEpsilonWallCondition2D);
}
}  // namespace Kratos.
