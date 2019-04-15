#include "custom_processes/define_wake_2d_process.h"

namespace Kratos
{
    void DefineWake2DProcess::Execute()
    {
        KRATOS_TRY;
        SaveTrailingEdgeNodecpp();
        MarkWakeElementscpp();


        KRATOS_CATCH("");
    }

    void DefineWake2DProcess::SaveTrailingEdgeNodecpp()
    {
        double max_x_coordinate = -1e30;
        auto trailing_edge_node = mrBodyModelPart.NodesBegin();
        for(auto it=mrBodyModelPart.NodesBegin(); it!=mrBodyModelPart.NodesEnd(); ++it)
        {
            if (it->X()>max_x_coordinate)
            {
                max_x_coordinate = it->X();
                trailing_edge_node = it;
            }
        }
        trailing_edge_node->SetValue(TRAILING_EDGE,true);
    }

    //void DefineWake2DProcess::MarkTrailingEdgeElementscpp(ModelPart& trailing_edge_model_part, int j)
    template <typename GE, typename IT> void DefineWake2DProcess::MarkTrailingEdgeElementscpp(ModelPart& trailing_edge_model_part, GE geom, IT it_elem)
    {
        // auto it_elem = mrFluidModelPart.ElementsBegin() + j;
        // auto geom = it_elem->GetGeometry();
        for (unsigned int k = 0; k < geom.size(); ++k)
        {
            if (geom[k].GetValue(TRAILING_EDGE))
            {
                it_elem->SetValue(TRAILING_EDGE, true);
                std::cout << "it_elem->" << it_elem->Id() << std::endl;
                //trailing_edge_model_part.Elements()(it_elem);
                auto elem = mrFluidModelPart.Elements()(it_elem->Id());
                trailing_edge_model_part.AddElement(elem);
                break;
            }
        }
    }

    template <typename GE, typename IT> void DefineWake2DProcess::SelectPotentiallyWakeElementscpp(GE geom, IT it_elem)
    {
        std::cout << "Hi" << std::endl;
    }
}