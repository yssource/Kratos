// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================
//

#if !defined(KRATOS_CUSTOM_WHOLE_CUTTING_PROCESS_H_INCLUDED )
#define  KRATOS_CUSTOM_WHOLE_CUTTING_PROCESS_H_INCLUDED


// System includes
// Please put system includes in custom_whole_cutting_process.h

// External includes
// Please put external includes in custom_whole_cutting_process.h

// Project includes

#include "custom_processes/custom_whole_cutting_process.h"

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include "math.h"

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
// Application includes
#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h" // To find node neighbours using elements
#include "processes/find_conditions_neighbours_process.h" // To find node neighbours using conditions
#include "processes/node_erase_process.h" // To delete empty nodes
#include "utilities/normal_calculation_utils.h" // To calculate element's normal
#include "geometries/triangle_3d_3.h" // Skin face geometry template


namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.

class CustomWholeCuttingProcess{
public:
	// Needed structures for the ExtractSurfaceMesh operation
	struct KeyComparor
	{
		bool operator()(const vector<unsigned int>& lhs, const vector<unsigned int>& rhs) const
		{
			if(lhs.size() != rhs.size())
				return false;

			for(unsigned int i=0; i<lhs.size(); i++)
			{
				if(lhs[i] != rhs[i]) return false;
			}

			return true;
		}
	};

	struct KeyHasher
	{
		std::size_t operator()(const vector<int>& k) const
		{
			return boost::hash_range(k.begin(), k.end());
		}
	};


	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomWholeCuttingProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomWholeCuttingProcess);

	///@}
	///@name Life Cycle
	///@{

	CustomWholeCuttingProcess() {
	}

	/// Destructor.
	virtual ~CustomWholeCuttingProcess() {
	}

	///@}
	///@name Operators
	///@{

	void operator()() {
		Execute();
	}

	///@}
	///@name Operations
	///@{

	/// For CHIMERA boundary condition purposes: Extracts a volume mesh with a certain threshold value
	void ExtractSurfaceMeshAtDistance(ModelPart& rModelPart, ModelPart& rExtractedSurfaceModelPart, double distance)
	{

		KRATOS_TRY;

		std::cout<<"\n::[Volume Mesh Extraction]::"<<std::endl;
		ModelPart rExtractedVolumeModelPart;

		// Initializing mesh nodes
		rExtractedVolumeModelPart.Nodes() = rModelPart.Nodes();

		// Extracting mesh elements which are only above the threshold value
		std::cout<<"  Extracting elements with in a distance of > " << fabs(distance) <<std::endl;
		Element::Pointer pElem;
		for(ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			double elementDistance = 0.0;
			unsigned int j = 0;
			for (j = 0 ; j < it->GetGeometry().PointsNumber(); j++){
				elementDistance += it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);
			}
			// Element distance is the average of the nodal distances.
			elementDistance = (elementDistance/j);

			if(elementDistance < distance)
			{
				pElem = Element::Pointer(new Element(*it));
				rExtractedVolumeModelPart.Elements().push_back(pElem);
			}
		}

		std::cout<<"  Successful extraction of the Volume !! "<<rExtractedVolumeModelPart.GetMesh()<<"\b"<<std::endl;
		ExtractSurfaceMesh(rExtractedVolumeModelPart, rExtractedSurfaceModelPart);

		KRATOS_CATCH("");

	}


	/// For CHIMERA boundary condition purposes: Extracts a volume mesh with a certain threshold value
	void ExtractVolumeMeshBetweenLimits(ModelPart& rModelPart,ModelPart& rExtractedVolumeModelPart, double lLimit, double uLimit)
	{

		KRATOS_TRY;

		std::cout<<"\n::[Volume Mesh Extraction]::"<<std::endl;

		// Initializing mesh nodes
		rExtractedVolumeModelPart.Nodes() = rModelPart.Nodes();

		// Extracting mesh elements which are only above the threshold value
		std::cout<<"  Extracting elements between " << lLimit <<" and "<<uLimit<<std::endl;
		Element::Pointer pElem;
		for(ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			double elementDistance = 0.0;
			int numPointsInside = 0;
			unsigned int j = 0;
			for (j = 0 ; j < it->GetGeometry().PointsNumber(); j++){
				elementDistance = it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);
				if(elementDistance > lLimit && elementDistance < uLimit )
				{
				  numPointsInside++;
				}
			}

			if(numPointsInside > 0)
			{
				pElem = Element::Pointer(new Element(*it));
				rExtractedVolumeModelPart.Elements().push_back(pElem);
			}
		}

		std::cout<<" ########  Successful extraction of the Volume !! "<<rExtractedVolumeModelPart.GetMesh()<<"\b"<<std::endl;
		KRATOS_CATCH("");
	}



	/// For Topology Optimization purposes: Extracts a surface mesh from a provided volume mesh
	void ExtractSurfaceMesh(ModelPart& rExtractedVolumeModelPart, ModelPart& rExtractedSurfaceModelPart)
	{
		KRATOS_TRY;

		std::cout<<"::[Surface Mesh Extraction]::"<<std::endl;

		// Some type-definitions
		typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasher, KeyComparor > hashmap;
		typedef boost::unordered_map<vector<unsigned int>, vector<unsigned int>, KeyHasher, KeyComparor > hashmap_vec;

		// Create map to ask for number of faces for the given set of node ids representing on face in the model part
		hashmap n_faces_map;

		// Fill map that counts number of faces for given set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for(unsigned int face=0; face<faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<unsigned int> ids(faces[face].size());

				// Store node ids
				for(unsigned int i=0; i<faces[face].size(); i++)
					ids[i] = faces[face][i].Id();

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				// Fill the map
				n_faces_map[ids] += 1;
			}
		}

		// Create a map to get nodes of skin face in original order for given set of node ids representing that face
		// The given set of node ids may have a different node order
		hashmap_vec ordered_skin_face_nodes_map;

		// Fill map that gives original node order for set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for(unsigned int face=0; face<faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<unsigned int> ids(faces[face].size());
				vector<unsigned int> unsorted_ids(faces[face].size());

				// Store node ids
				for(unsigned int i=0; i<faces[face].size(); i++)
				{
					ids[i] = faces[face][i].Id();
					unsorted_ids[i] = faces[face][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				if(n_faces_map[ids] == 1)
					ordered_skin_face_nodes_map[ids] = unsorted_ids;
			}
		}
		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		unsigned int id_condition = 1;
		rExtractedSurfaceModelPart.Nodes() = rExtractedVolumeModelPart.Nodes();

		// Add skin faces as triangles to skin-model-part (loop over all node sets)
		std::cout<<"  Extracting surface mesh and computing normals" <<std::endl;
		for(typename hashmap::const_iterator it=n_faces_map.begin(); it!=n_faces_map.end(); it++)
		{
			// If given node set represents face that is not overlapping with a face of another element, add it as skin element
			if(it->second == 1)
			{
				// If skin face is a triangle store triangle in with its original orientation in new skin model part
				if(it->first.size()==3)
				{
					// Getting original order is important to properly reproduce skin face including its normal orientation
					vector<unsigned int> original_nodes_order = ordered_skin_face_nodes_map[it->first];
					Node < 3 >::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node < 3 >::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node < 3 >::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const& rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Skin faces are added as conditions
					Triangle3D3< Node<3> > triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);
				}
				// If skin face is a quadrilateral then divide in two triangles and store them with their original orientation in new skin model part
				if(it->first.size()==4)
				{
					// Getting original order is important to properly reproduce skin including its normal orientation
					vector<unsigned int> original_nodes_order = ordered_skin_face_nodes_map[it->first];

					Node < 3 >::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node < 3 >::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node < 3 >::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					Node < 3 >::Pointer pnode4 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[3]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const& rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Add triangle one as condition
					Triangle3D3< Node<3> > triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);

					// Add triangle two as condition
					Triangle3D3< Node<3> > triangle2(pnode1, pnode3, pnode4);
					Condition::Pointer p_condition2 = rReferenceTriangleCondition.Create(id_condition++, triangle2, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition2);
				}
			}
		}
		std::cout<<"Successful extraction of the Surface "<<rExtractedSurfaceModelPart.GetMesh()<<std::endl;

		KRATOS_CATCH("");
	}






	virtual void Execute() {
	}

	virtual void Clear() {
	}

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const {
		return "CustomWholeCuttingProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << "CustomWholeCuttingProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
	}







	///@}
	///@name Friends
	///@{

	///@}

protected:

	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:

	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	CustomWholeCuttingProcess& operator=(CustomWholeCuttingProcess const& rOther);

	/// Copy constructor.
	//CustomWholeCuttingProcess(CustomWholeCuttingProcess const& rOther);

	///@}

}; // Class CustomWholeCuttingProcess

}  // namespace Kratos.

#endif // KRATOS_CUSTOM_WHOLE_CUTTING_PROCESS_H_INCLUDED  defined
