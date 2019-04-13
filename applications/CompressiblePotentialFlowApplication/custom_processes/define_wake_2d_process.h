#ifndef KRATOS_DEFINE_WAKE_2D_PROCESS_H
#define KRATOS_DEFINE_WAKE_2D_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class DefineWake2DProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DefineWake2DProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;
    ///@}
    ///@name Life Cycle
    ///@{

    DefineWake2DProcess(ModelPart& rModelPart
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~DefineWake2DProcess() override {}



    ///@}
    ///@name Operators
    ///@{
    void PrintMyInfo2();
    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }



    // void PrintMyInfo2()
    // {
    //     KRATOS_TRY
    //     std::cout<< "Print My Info !!!!!!!!" << std::endl;
    //     KRATOS_CATCH("");
    // }

    ///@}
    ///@name Operations
    ///@{

    //void SaveTrailingEdgeNodecpp();
    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;
        //SaveTrailingEdgeNodecpp();
        PrintMyInfo2();
        double max_x_coordinate = -1e30;
        auto trailing_edge_node = mrModelPart.NodesBegin();
        for(auto it=mrModelPart.NodesBegin(); it!=mrModelPart.NodesEnd(); ++it)
        //for(auto it_elem = mrModelPart.ElementsBegin(); it_elem!=mrModelPart.ElementsEnd(); ++it_elem) // Loop the elements
        {
            if (it->X()>max_x_coordinate)
            {
                max_x_coordinate = it->X();
                trailing_edge_node = it;
            }
        }
        std::cout<< "max_x_coordinate cpp" <<  max_x_coordinate << std::endl;
        std::cout<< "trailing_edge_node" <<  trailing_edge_node->Id() << std::endl;
        int temp;
        std::cin>>temp;
        //it->Set(TRAILING_EDGE, max_x_coordinate); // how to save it?

        KRATOS_CATCH("");
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
        return "DefineWake2DProcess";
    }


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DefineWake2DProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


private:

    ModelPart& mrModelPart;
    Flags mrOptions;


    /// Assignment operator.
    DefineWake2DProcess& operator=(DefineWake2DProcess const& rOther);

    /// Copy constructor.
    DefineWake2DProcess(DefineWake2DProcess const& rOther);
    //void PrintMyInfo2();

}; // Class Process



/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DefineWake2DProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DefineWake2DProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // KRATOS_DEFINE_WAKE_2D_PROCESS_H