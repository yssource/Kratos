//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


// System includes
#include <chrono>

// External includes


// Project includes
#include "processes/calculate_distance_to_terrain_process.h"


namespace Kratos
{

CalculateDistanceToTerrainProcess::CalculateDistanceToTerrainProcess(
    ModelPart &rVolumePart, ModelPart &rTerrainPart, Parameters TheParameters)
    : CalculateDistanceToSkinProcess(rVolumePart, rTerrainPart, TheParameters) {
}

CalculateDistanceToTerrainProcess::~CalculateDistanceToTerrainProcess() {}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	double CalculateDistanceToTerrainProcess::DistancePositionInSpace(double* coords)
	{

		typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

        intersections_container_type intersections;

        const int dimension = 3;

        double distances[3] = {1.0, 1.0, 1.0};
        bool ray_is_valid[3]={true,true,true};
        
        std::size_t z_direction = 2;

        // Creating the ray
        double ray[3] = {coords[0], coords[1], coords[2]};

        OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess::CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();
        pOctree->NormalizeCoordinates(ray);
        ray[z_direction] = 0; // starting from the lower extreme

        this->GetRayIntersections(ray, z_direction, intersections);

        if(intersections.empty())
            return 1.00; // I'm not sure if returning max double would cause a problem. Pooyan.

        auto last_hit = intersections.back().first;
        return coords[z_direction] - last_hit;

	}

	/// Turn back information as a string.
	std::string CalculateDistanceToTerrainProcess::Info() const
	{
		return "CalculateDistanceToTerrainProcess";
	}

	/// Print information about this object.
	void CalculateDistanceToTerrainProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void CalculateDistanceToTerrainProcess::PrintData(std::ostream& rOStream) const
	{
	}



}  // namespace Kratos.
