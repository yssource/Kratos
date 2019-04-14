#ifndef KRATOS_DEFINE_WAKE_2D_PROCESS_H
#define KRATOS_DEFINE_WAKE_2D_PROCESS_H

// #include "includes/define.h"
// #include "includes/model_part.h"
// #include "includes/kratos_flags.h"
// #include "processes/process.h"
// #include "geometries/geometry.h"
// #include "utilities/geometry_utilities.h"
#include "compressible_potential_flow_application_variables.h"
// #include "utilities/math_utils.h"
// #include "includes/kratos_parameters.h"

// #include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
// #include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
// #include <utility>

#include <string>
#include <iostream>
#include <sstream>

#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

//class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) DefineWake2DProcess: public Process
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

    // DefineWake2DProcess(
    //     ModelPart& rModelPart,
    //     Parameters Settings
    //     ) : Process(Flags()) ,
    //         mrModelPart(rModelPart),
    //         mSettings( Settings)
    // {
    //     KRATOS_TRY
    //     Parameters default_parameters( R"(
    //     {
    //         "model_part_name": "",
    //         "wake_direction": [1.0,0.0,0.0],
    //         "epsilon": 1e-9
    //     }  )" );

    //     Settings.ValidateAndAssignDefaults(default_parameters);
    //     KRATOS_CATCH("")
    // }

    DefineWake2DProcess(
        ModelPart& rModelPart
        ) : Process(Flags()) ,
            mrModelPart(rModelPart)
    {
    }

    /// Assignment operator.
    DefineWake2DProcess& operator=(DefineWake2DProcess const& rOther);
    /// Destructor.
    ~DefineWake2DProcess() override = default;
    /// Copy constructor.
    DefineWake2DProcess(DefineWake2DProcess const& rOther);




    ///@}
    ///@name Operators
    ///@{
    //void PrintMyInfo2();
    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    //void SaveTrailingEdgeNodecpp();
    //void MarkWakeElementscpp();
    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;
        SaveTrailingEdgeNodecpp();
        //MarkWakeElementscpp();


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
    //Parameters mSettings;   /// The settings of the problem (names of the conditions and elements)
    Flags mrOptions;

    void MarkWakeElementscpp()
    {
        //ModelPart* ptrailing_edge_model_part = &mrModelPart.CreateSubModelPart("fluid_model_part");
        //trailing_edge_model_part->CreateSubModelPart("trailing_edge_model_part");
        //ModelPart trailing_edge_model_part = mrModelPart.CreateSubModelPart("trailing_edge_model_part");
        for(auto it_elem = mrModelPart.ElementsBegin(); it_elem!=mrModelPart.ElementsEnd(); ++it_elem) // Loop the elements
        {
            auto geom = it_elem->GetGeometry();

            for(unsigned int i=0; i<geom.size(); ++i)
            {
                 if (geom[i].GetValue(TRAILING_EDGE))
                 {
                     it_elem->SetValue(TRAILING_EDGE,true);
                     std::cout<<"it_elem->"<<it_elem->Id()<<std::endl;
                     break;
                 }
            }
            //trailing_edge_model_part.Elements()(it_elem);// how to add the elements to the model part

            //MarkTrailingEdgeElementscpp(it_elem);
        }
    }

    void SaveTrailingEdgeNodecpp()
    {
        double max_x_coordinate = -1e30;
        auto trailing_edge_node = mrModelPart.NodesBegin();
        //trailing_edge_node->Set(TRAILING_EDGE,1500); // how to save it?
        for(auto it=mrModelPart.NodesBegin(); it!=mrModelPart.NodesEnd(); ++it)
        {
            if (it->X()>max_x_coordinate)
            {
                max_x_coordinate = it->X();
                trailing_edge_node = it;
                //it->Set(TRAILING_EDGE,true)
            }
        }
        trailing_edge_node->SetValue(TRAILING_EDGE,true);
        std::cout<< "max_x_coordinate cpp: " <<  max_x_coordinate << std::endl;
        std::cout<< "trailing_edge_node cpp: " <<  trailing_edge_node->Id() << std::endl;
        //std::cout<< "TRAILING_EDGE cpp: " <<  TRAILING_EDGE << std::endl;
    }

    // void MarkTrailingEdgeElementscpp(auto it_elem)
    // {
    //     std::cout<< "elem nodes " <<it_elem.node<<std::endl;
    // }




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