//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ERODIBLE_BED_APPLICATION_H_INCLUDED )
#define  KRATOS_ERODIBLE_BED_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
#include "custom_elements/erodible_bed_1d.h" //including the file for the element
#include "custom_elements/erodible_bed_2d.h" //including the file for the element
#include "custom_utilities/bed_particle.h" //hahaha

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	KRATOS_DEFINE_VARIABLE(double, SEDIMENT_VELOCITY_OVER_ELEM_SIZE)
	KRATOS_DEFINE_VARIABLE(double, MEAN_SIZE)
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_BED_PARTICLES)
	KRATOS_DEFINE_VARIABLE(int, BED_PARTICLE_POINTERS_OFFSET)	
	KRATOS_DEFINE_VARIABLE(double, HEIGHT)
	KRATOS_DEFINE_VARIABLE(double, PROJECTED_HEIGHT)
	KRATOS_DEFINE_VARIABLE(double, DELTA_HEIGHT)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SEDIMENT_VELOCITY)


	typedef PointerVector< Bed_Particle, Bed_Particle*, std::vector<Bed_Particle*> >BedParticlePointerVector;
	KRATOS_DEFINE_VARIABLE( BedParticlePointerVector , BED_PARTICLE_POINTERS)	


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
	/** Detail class definition.
	*/
	class KratosErodibleBedApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosErodibleBedApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosErodibleBedApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosErodibleBedApplication();

		/// Destructor.
		virtual ~KratosErodibleBedApplication(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



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
			return "KratosErodibleBedApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in my application");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
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



		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
		///@name Member Variables 
		///@{ 
// 		const Elem2D   mElem2D; 
 		const ErodibleBed1D   mErodibleBed1D; 
 		const ErodibleBed2D   mErodibleBed2D; 

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
		KratosErodibleBedApplication& operator=(KratosErodibleBedApplication const& rOther);

		/// Copy constructor.
		KratosErodibleBedApplication(KratosErodibleBedApplication const& rOther);


		///@}    

	}; // Class KratosErodibleBedApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_ERODIBLE_BED_APPLICATION_H_INCLUDED  defined 


