#include "custom_processes/define_wake_2d_process.h"

namespace Kratos
{
    void DefineWake2DProcess::Execute()
    {
        KRATOS_TRY;
        auto trailing_edge = mrBodyModelPart.NodesBegin();
        SaveTrailingEdgeNodecpp(trailing_edge);
        MarkWakeElementscpp(trailing_edge);


        KRATOS_CATCH("");
    }

    //void DefineWake2DProcess::SaveTrailingEdgeNodecpp()
    template <typename TE> void DefineWake2DProcess::SaveTrailingEdgeNodecpp(TE &trailing_edge_node)
    {
        double max_x_coordinate = -1e30;
        //auto trailing_edge_node = mrBodyModelPart.NodesBegin();
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

    template <typename GE, typename IT, typename TE> void DefineWake2DProcess::SelectPotentiallyWakeElementscpp(GE geom, IT it_elem, TE trailing_edge)
    {
        double x_distance_to_te = it_elem->GetGeometry().Center().X() - trailing_edge->X();
        double y_distance_to_te = it_elem->GetGeometry().Center().Y() - trailing_edge->Y();
    }

    template <typename TE> void DefineWake2DProcess::MarkWakeElementscpp(TE trailing_edge)
    {
        //ModelPart* ptrailing_edge_model_part = &mrModelPart.CreateSubModelPart("fluid_model_part");
        //trailing_edge_model_part->CreateSubModelPart("trailing_edge_model_part");
        ModelPart& trailing_edge_model_part = mrFluidModelPart.CreateSubModelPart("trailing_edge_model_part2");

        for(int i = 0; i< static_cast<int>(mrFluidModelPart.Elements().size()); ++i) // Loop the elements
        {
            auto it_elem = mrFluidModelPart.ElementsBegin() + i;
            auto geom = it_elem->GetGeometry();
            MarkTrailingEdgeElementscpp(trailing_edge_model_part, geom, it_elem);
            SelectPotentiallyWakeElementscpp(geom, it_elem, trailing_edge);
        }
    }
}