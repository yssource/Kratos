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

#if !defined(KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED)
#define KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED

// System includes

#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "utilities/binbased_fast_point_locator.h"

// Project includes

#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

namespace Kratos
{

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
template <unsigned int TDim>
class CustomExtractVariablesProcess
{
  public:
	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomExtractVariablesProcess);

	///@}
	///@name Life Cycle
	///@{
	typedef array_1d<double, TDim + 1> VariableArrayType;
	typedef Geometry<Node<3>>::CoordinatesArrayType CoordinatesArrayType;
	typedef typename BinBasedFastPointLocator<TDim>::Pointer BinBasedPointLocatorPointerType;

	CustomExtractVariablesProcess()
	{
		//this->pBinLocator = NULL;
	}

	/// Destructor.
	virtual ~CustomExtractVariablesProcess()
	{
	}

	///@}
	///@name Operators
	///@{

	void operator()()
	{
		Execute();
	}

	///@}
	///@name Operations
	///@{

	virtual void Execute()
	{
	}

	virtual void Clear()
	{
	}
	// Extracts the velocity components and pressure at any point

	void ExtractVariables(ModelPart &fromVolumeModelPart,double X,double Y,double Z, double* rvariables)

	{
		array_1d<double,3> rcoord ;
		rcoord[0] = X;
		rcoord[1] = Y;
		rcoord[2] = Z;
		
		BinBasedPointLocatorPointerType pBinLocator = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(fromVolumeModelPart));

		array_1d<double, TDim + 1> N;
		const int max_results = 10000;
		typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

		typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		Element::Pointer pElement;
		bool is_found = false;
		is_found = pBinLocator->FindPointOnMesh(rcoord, N, pElement, result_begin, max_results);

		if (is_found == true)
		{

			Geometry<Node<3>> &geom = pElement->GetGeometry();

			if (TDim == 2)
			{
				for (int i = 0; i < TDim + 1; i++)
				{
					
					rvariables[0] += N[i] * geom[i].FastGetSolutionStepValue(VELOCITY_X);
					rvariables[1] += N[i] * geom[i].FastGetSolutionStepValue(VELOCITY_Y);
					rvariables[2] += N[i] * geom[i].FastGetSolutionStepValue(PRESSURE);
				}
				if (TDim == 3)
				{
					for (int i = 0; i < TDim + 1; i++)
					{

						rvariables[0] += N[i] * geom[i].FastGetSolutionStepValue(VELOCITY_X);
						rvariables[1] += N[i] * geom[i].FastGetSolutionStepValue(VELOCITY_Y);
						rvariables[2] += N[i] * geom[i].FastGetSolutionStepValue(VELOCITY_Z);
						rvariables[3] += N[i] * geom[i].FastGetSolutionStepValue(PRESSURE);
					}
				}
			}
		}

		
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
	virtual std::string Info() const
	{
		return "CustomExtractVariablesProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "CustomExtractVariablesProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
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

	//ModelPart &mrBackGroundModelPart;
	//ModelPart &mrPatchSurfaceModelPart;
	//BinBasedFastPointLocator<TDim> *pBinLocator; // Template argument 3 stands for 3D case
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
	CustomExtractVariablesProcess &operator=(CustomExtractVariablesProcess const &rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}

}; // Class CustomExtractVariablesProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
