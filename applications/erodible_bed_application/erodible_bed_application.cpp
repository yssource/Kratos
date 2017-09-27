//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "erodible_bed_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
	KRATOS_CREATE_VARIABLE(double, SEDIMENT_VELOCITY_OVER_ELEM_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_SIZE)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_BED_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, BED_PARTICLE_POINTERS_OFFSET)	
	KRATOS_CREATE_VARIABLE(double, HEIGHT)
	KRATOS_CREATE_VARIABLE(double, PROJECTED_HEIGHT)
	KRATOS_CREATE_VARIABLE(double, DELTA_HEIGHT)//
	KRATOS_CREATE_VARIABLE( BedParticlePointerVector , BED_PARTICLE_POINTERS)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SEDIMENT_VELOCITY)



 	KratosErodibleBedApplication::KratosErodibleBedApplication(): //constructor  do not forget to add the ":" 
		mErodibleBed1D    ( 0, Element::GeometryType::Pointer( new Line2D<Node<3> >(  Element::GeometryType::PointsArrayType (2 ) ) ) ),
		mErodibleBed2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) )
		{}
 	
 	void KratosErodibleBedApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosErodibleBedApplication... " << std::endl;
	 
		KRATOS_REGISTER_VARIABLE(SEDIMENT_VELOCITY_OVER_ELEM_SIZE)
		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
		KRATOS_REGISTER_VARIABLE(NUMBER_OF_BED_PARTICLES)
		KRATOS_REGISTER_VARIABLE(BED_PARTICLE_POINTERS_OFFSET)	
		KRATOS_REGISTER_VARIABLE(HEIGHT)
		KRATOS_REGISTER_VARIABLE(PROJECTED_HEIGHT)
		KRATOS_REGISTER_VARIABLE(DELTA_HEIGHT)//
		KRATOS_REGISTER_VARIABLE(BED_PARTICLE_POINTERS)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SEDIMENT_VELOCITY)


		KRATOS_REGISTER_ELEMENT("ErodibleBed1D", mErodibleBed1D);  //and here is our element
		KRATOS_REGISTER_ELEMENT("ErodibleBed2D", mErodibleBed2D);  //and here is our element

 
 	}

}  // namespace Kratos.


