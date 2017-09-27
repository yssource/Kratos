//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: ---$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_DEM_utils_application.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "geometries/sphere_3d_1.h"
//#include "../DEM_application/DEM_application.h"

namespace Kratos
{
        
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY) 

  
KratosCustomDEMUtilsApplication::KratosCustomDEMUtilsApplication()//:
  //mSwimmingNanoParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
  //mSwimmingAnalyticParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1))))
{}

void KratosCustomDEMUtilsApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosCustomDEMUtilsApplication... " << std::endl;
                
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)   

  /* Define In Global variables.cpp */

  //KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition2D",  mComputeLaplacianSimplexCondition2D)
  //KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition3D", mComputeLaplacianSimplexCondition3D)
  Serializer::Register( "VariablesList", mVariablesList );
 }

}  // namespace Kratos.


