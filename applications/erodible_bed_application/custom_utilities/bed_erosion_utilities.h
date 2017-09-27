

#if !defined(KRATOS_ADD_BED_EROSION_UTILITY_INCLUDED )
#define KRATOS_ADD_BED_EROSION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "erodible_bed_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class BedErosionUtility
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(BedErosionUtility);

		BedErosionUtility(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~BedErosionUtility()
		{}

		
		void CalculateSedimentVelocity() 
		{
			KRATOS_TRY
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				inode->FastGetSolutionStepValue(SEDIMENT_VELOCITY) = 0.1 * inode->FastGetSolutionStepValue(VELOCITY);
			}
			
			KRATOS_CATCH("")
		} 
		
		void MoveBedMesh(double depth)
		{
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				(inode)->Z() = (inode)->Z0() - depth + inode->GetSolutionStepValue(HEIGHT);
			}
		}

		void MoveMeshUsingUserVariable(ModelPart::NodesContainerType& rNodes,const Variable< array_1d<double, 3 > >& MyVariable)
		{
			for(ModelPart::NodesContainerType::iterator inode = rNodes.begin(); 
				inode!=rNodes.end(); inode++)

			{
				const array_1d<double,3> displacement =  inode->FastGetSolutionStepValue(MyVariable);
				(inode)->X() = (inode)->X0() + displacement[0];
				(inode)->Y() = (inode)->Y0() + displacement[1];
				(inode)->Z() = (inode)->Z0() + displacement[2];
			}
		}

	protected:


	private:
		ModelPart& mr_model_part;

	};
	
	
}  // namespace Kratos.

#endif // KRATOS_ADD_BED_EROSION__UTILITY_INCLUDED  defined 


