#ifndef POSTPROCESS_UTILITIES_H
#define POSTPROCESS_UTILITIES_H

#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/variables.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{

class PostProcessUtilities
{

public:
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

    /// Default constructor.

    PostProcessUtilities(){};

    /// Destructor.

    virtual ~PostProcessUtilities(){};

    void RebuildPostProcessModelPart(ModelPart &r_post_model_part, ModelPart &r_main_model_part)
    {
        r_post_model_part.Elements().clear();
        r_post_model_part.Nodes().clear();
        for (size_t i = 0; i < r_main_model_part.NumberOfNodes(); i++)
        {
            auto node = r_main_model_part.NodesBegin() + i;
            r_post_model_part.AddNode(*(node.base()));
        }

        PointerVector<Element> elements_to_be_added_first;
        PointerVector<Element> elements_to_be_added_later;

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_main_model_part.SubModelPartsBegin();
             sub_model_part != r_main_model_part.SubModelPartsEnd(); ++sub_model_part)
        {

            ModelPart &smp_k = *sub_model_part;
            for (size_t i = 0; i < smp_k.NumberOfElements(); i++)
            {
                auto elem = smp_k.ElementsBegin() + i;
                if (r_main_model_part.GetMesh(0).HasElement(elem->Id()))
                {
                    if (r_main_model_part.GetElement(elem->Id()).GetGeometry().GetGeometryType() == elem->GetGeometry().GetGeometryType())
                    {
                        elements_to_be_added_first.push_back(*(elem.base()));
                    }
                    else
                    {
                        elements_to_be_added_later.push_back(*(elem.base()));
                    }
                }
                else
                {
                    elements_to_be_added_later.push_back(*(elem.base()));
                }
            }
        }
        r_post_model_part.AddElements(elements_to_be_added_first.begin(), elements_to_be_added_first.end());

        int max_id = 1;
        for (size_t i = 0; i < r_post_model_part.NumberOfElements(); i++)
        {
            auto elem = r_post_model_part.ElementsBegin() + i;
            if ((int)elem->Id() > max_id)
                max_id = elem->Id();
        }

        for (auto &elem : elements_to_be_added_later)
        {
            max_id += 1;
            elem.SetId(max_id);
            elem.Set(ACTIVE, true);
        }
        r_post_model_part.AddElements(elements_to_be_added_later.begin(), elements_to_be_added_later.end());
    }

protected:
}; // Class PostProcessUtilities

} // namespace Kratos.

#endif // POSTPROCESS_UTILITIES_H
