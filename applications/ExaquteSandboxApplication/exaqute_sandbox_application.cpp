//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes


// External includes


// Project includes
#include "exaqute_sandbox_application.h"
#include "exaqute_sandbox_application_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"


namespace Kratos {

KratosExaquteSandboxApplication::KratosExaquteSandboxApplication():
    KratosApplication("ExaquteSandboxApplication"),
    mFractionalStepSemiExplicit2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFractionalStepSemiExplicit3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
    {}

void KratosExaquteSandboxApplication::Register()
{

  // calling base class register to register Kratos components
  KratosApplication::Register();
  KRATOS_INFO("") << "Initializing KratosExaquteSandboxApplication..." << std::endl;
  // varaibles
  KRATOS_REGISTER_VARIABLE( DIVERGENCE_WEIGHTED )
  KRATOS_REGISTER_VARIABLE( VELOCITY_H1_SEMINORM )
  // elements
  KRATOS_REGISTER_ELEMENT("FractionalStepSemiExplicitElement2D3N", mFractionalStepSemiExplicit2D3N);
  KRATOS_REGISTER_ELEMENT("FractionalStepSemiExplicitElement3D4N", mFractionalStepSemiExplicit3D4N);

}
}  // namespace Kratos.
