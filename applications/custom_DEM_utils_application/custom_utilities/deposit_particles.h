#if !defined(KRATOS_ADD_DEPOSIT_PARTICLES_UTILITY_INCLUDED )
#define KRATOS_ADD_DEPOSIT_PARTICLES_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "includes/define.h"
#include "custom_DEM_utils_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "../DEM_application/DEM_application_variables.h"
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
	class DepositParticlesUtility
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(DepositParticlesUtility);

		DepositParticlesUtility(ModelPart& walls_model_part)
			: mr_walls_model_part(walls_model_part)
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;


			if(mr_walls_model_part.NodesBegin()->SolutionStepsDataHas(THICKNESS) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing THICKNESS variable on wall model part (could not find it when initializing DepositParticlesUtility","");
			if(mr_walls_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on wall model part (could not find it when initializing DepositParticlesUtility","");


			KRATOS_CATCH("")	
		}
		

		~DepositParticlesUtility()
		{}

		
		void IncreaseNodalThicknessUsingParticles(bool first_dem_iter) 
		{
			KRATOS_TRY
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			if(first_dem_iter)
				for(ModelPart::NodesContainerType::iterator inode = mr_walls_model_part.NodesBegin(); 
					inode!=mr_walls_model_part.NodesEnd(); inode++)
				{
					inode->FastGetSolutionStepValue(NODAL_AREA) =0.0;
					inode->FastGetSolutionStepValue(THICKNESS) =0.0;
				}

			for(ModelPart::ConditionsContainerType::iterator icond = mr_walls_model_part.ConditionsBegin(); 
				icond!=mr_walls_model_part.ConditionsEnd(); icond++)
			{	
				const int number_of_nodes=3;				
				//before proceeding, we make sure we are a "wall"(bed) with friciton (that is, not a slip wall) 
				bool bed=true;
				if (icond->GetProperties()[WALL_FRICTION]<1.0e-3)
					bed=false;
				if(bed)
				{
					Geometry< Node<3> >& geom = icond->GetGeometry();
					
					DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*icond));

					std::vector<SphericParticle*>& rNeighbours = p_wall->mNeighbourSphericParticles;
				
					for (unsigned int i=0; i<rNeighbours.size(); i++)
					{
						
						const array_1d<double,3> particle_coord = rNeighbours[i]->GetGeometry()[0].Coordinates();
							
						array_1d<double,number_of_nodes> N;
						double Area=0.0;
						bool is_found = CalculatePosition(geom,particle_coord[0],particle_coord[1],N,Area);
						if(is_found)
						{
							rNeighbours[i]->Set(TO_ERASE);
							//COMPUTE HERE MASS; VOLUME; WHATEVER!
							const double this_particle_volume = 0.000000001;
							const double this_particle_thickness_contribution =  this_particle_volume/Area ;
							for (unsigned int k=0; k< number_of_nodes; k++)
							{

								geom[k].FastGetSolutionStepValue(THICKNESS) += this_particle_thickness_contribution*N[k] * Area;
							}
						}
					
					}//Loop condition neighbours (spheres)
					
					//assembling nodal area
					if (first_dem_iter)
					{
						array_1d<double,number_of_nodes> dummy_N;
						double Area=0.0;
						CalculatePosition(geom,0.0,0.0,dummy_N,Area);

						for (unsigned int k=0; k< number_of_nodes; k++)
						{
							geom[k].FastGetSolutionStepValue(NODAL_AREA) +=  1.0/3.0 * Area;
						}
					}
				}
			} //loop on conditions
			KRATOS_CATCH("")

		}
		
		void NormalizeThickness()
		{
			KRATOS_TRY

			//finally we can divide by the nodal area.
			for(ModelPart::NodesContainerType::iterator inode = mr_walls_model_part.NodesBegin(); 
				inode!=mr_walls_model_part.NodesEnd(); inode++)
			{
				double& thickness = inode->FastGetSolutionStepValue(THICKNESS);
				const double nodal_area = inode->FastGetSolutionStepValue(NODAL_AREA);
				if(nodal_area>=1.0e-9)
					thickness /= inode->FastGetSolutionStepValue(NODAL_AREA);
								
			}
			KRATOS_CATCH("")

		}
		
	protected:

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc,
                array_1d<double, 3 > & N, double& total_triangle_area
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();

            total_triangle_area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (total_triangle_area == 0.0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area  = 1.0 / total_triangle_area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
			//KRATOS_WATCH(N);

            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }
        
        
        inline double CalculateVol(const double x0, const double y0,
                const double x1, const double y1,
                const double x2, const double y2
                )
        {
            return fabs(0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0)));
        }
        
	private:
		ModelPart& mr_walls_model_part;

	};
	
	
}  // namespace Kratos.

#endif // KRATOS_ADD_BED_EROSION__UTILITY_INCLUDED  defined 


