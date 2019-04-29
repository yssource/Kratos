//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_LOGARITHMIC_Y_PLUS_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LOGARITHMIC_Y_PLUS_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "rans_modelling_application_variables.h"

namespace Kratos
{
///@addtogroup RANSModellingApplication
///@{

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

/// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
/** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
    so that the fluid element can take natural convection into account.

    This process makes use of the following data:
    - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
    - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
    - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
    - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

    With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
    The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

    If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
    This is the usual value for perfect gases (if the temperature is given in Kelvin).
 */

class RansLogarithmicYPlusModelSensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::ElementType ElementType;

    typedef Geometry<NodeType> GeometryType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Pointer definition of RansLogarithmicYPlusModelSensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLogarithmicYPlusModelSensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansLogarithmicYPlusModelSensitivitiesProcess(ModelPart& rModelPart, Parameters& rParameters)
        : mrModelPart(rModelPart), mrParameters(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"      : 0,
            "step"            : 0,
            "constants": {
                "von_karman"  : 0.41,
                "beta"        : 5.2
            }
        })");

        mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mEchoLevel = mrParameters["echo_level"].GetInt();
        mStep = mrParameters["step"].GetInt();

        mVonKarman = mrParameters["constants"]["von_karman"].GetDouble();
        mBeta = mrParameters["constants"]["beta"].GetDouble();

        KRATOS_CATCH("");
    }
    /// Destructor.
    ~RansLogarithmicYPlusModelSensitivitiesProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS);
        KRATOS_CHECK_VARIABLE_KEY(RANS_Y_PLUS_SENSITIVITIES);

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
        }

        return 0;
    }

    void Execute() override
    {
        ModelPart::ElementsContainerType& r_elements = mrModelPart.Elements();

        const int number_of_elements = r_elements.size();

        const int domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

        const double inv_kappa = 1.0 / mVonKarman;

#pragma omp parallel for
        for (int i_element = 0; i_element < number_of_elements; ++i_element)
        {
            ElementType& r_element = *(r_elements.begin() + i_element);
            GeometryType& r_geometry = r_element.GetGeometry();
            const int number_of_nodes = r_geometry.PointsNumber();

            Matrix r_adjoint_y_plus_matrix(number_of_nodes, domain_size);

            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                const NodeType& r_node = r_geometry[i_node];
                const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
                const double wall_distance = r_node.FastGetSolutionStepValue(DISTANCE);
                const array_1d<double, 3> velocity =
                    r_node.FastGetSolutionStepValue(VELOCITY);
                const double velocity_magnitude = norm_2(velocity);
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

                double value = 0.0;

                if (y_plus > 11.06)
                {
                    value = (inv_kappa * std::log(y_plus) + mBeta);
                    value = value / (std::pow(value, 2) + velocity_magnitude * wall_distance *
                                                              inv_kappa / (nu * y_plus));
                }
                else
                {
                    value = 1.0 / (2.0 * y_plus);
                }
                for (int i_dim = 0; i_dim < domain_size; ++i_dim)
                {
                    r_adjoint_y_plus_matrix(i_node, i_dim) =
                        (wall_distance / nu) * value * velocity[i_dim] / velocity_magnitude;
                }
            }
            r_element.SetValue(RANS_Y_PLUS_SENSITIVITIES, r_adjoint_y_plus_matrix);
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
    std::string Info() const override
    {
        return std::string("RansLogarithmicYPlusModelSensitivitiesProcess");
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPart& mrModelPart;
    Parameters& mrParameters;

    unsigned int mEchoLevel;
    unsigned int mStep;

    double mVonKarman;
    double mBeta;

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
    RansLogarithmicYPlusModelSensitivitiesProcess& operator=(
        RansLogarithmicYPlusModelSensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansLogarithmicYPlusModelSensitivitiesProcess(RansLogarithmicYPlusModelSensitivitiesProcess const& rOther);

    ///@}

}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansLogarithmicYPlusModelSensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LOGARITHMIC_Y_PLUS_SENSITIVITIES_PROCESS_H_INCLUDED defined
