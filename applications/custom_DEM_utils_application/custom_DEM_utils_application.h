//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  ---$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_CUSTOM_DEM_UTILS_APPLICATION_H_INCLUDED )
#define  KRATOS_CUSTOM_DEM_UTILS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"  //TODO: must be removed eventually
#include "includes/legacy_structural_app_vars.h"  //TODO: must be removed eventually

#include "../DEM_application/custom_elements/spheric_particle.h"
#include "../DEM_application/DEM_application_variables.h"


namespace Kratos
{

    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)
  

class KratosCustomDEMUtilsApplication : public KratosApplication
{
public:

    /// Pointer definition of KratosCustomDEMUtilsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosCustomDEMUtilsApplication);

    /// Default constructor.
    KratosCustomDEMUtilsApplication();

    /// Destructor.
    virtual ~KratosCustomDEMUtilsApplication() {}


    virtual void Register();

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosCustomDEMUtilsApplication";
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
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

protected:

private:
    ///@name Static Member Variables
    ///@{
    //const SphericSwimmingParticle<NanoParticle> mSwimmingNanoParticle3D;
    //const SphericSwimmingParticle<AnalyticSphericParticle> mSwimmingAnalyticParticle3D;

    //const DEM_FEM_Particle mDEM_FEM_Particle2D;
    const VariablesList mVariablesList;

    ///@}
    ///@name Member Variables
    ///@{

    /// Assignment operator.
    KratosCustomDEMUtilsApplication& operator=(KratosCustomDEMUtilsApplication const& rOther);

    /// Copy constructor.
    KratosCustomDEMUtilsApplication(KratosCustomDEMUtilsApplication const& rOther);

}; // Class KratosCustomDEMUtilsApplication

}  // namespace Kratos.

#endif // KRATOS_CUSTOM_DEM_UTILS_APPLICATION_H_INCLUDED  defined 


