//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_POTENTIALFLOW_LIFT_OBJECTIVE_FUNCTION)
#define KRATOS_POTENTIALFLOW_LIFT_OBJECTIVE_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/objective_function.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for drag.
template <unsigned int TDim>
class PotentialFlowLiftObjectiveFunction : public ObjectiveFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PotentialFlowLiftObjectiveFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PotentialFlowLiftObjectiveFunction(Parameters& rParameters)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "objective_type": "drag",
            "structure_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "lift_direction": [0.0, 1.0, 0.0]
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mStructureModelPartName = rParameters["structure_model_part_name"].GetString();

        if (rParameters["lift_direction"].IsArray() == false ||
                rParameters["lift_direction"].size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "lift_direction vector is not a vector or does "
                               "not have size 3:",
                               rParameters.PrettyPrintJsonString())
        }

        for (unsigned int d = 0; d < TDim; ++d)
            mLiftDirection[d] = rParameters["lift_direction"][d].GetDouble();

        if (std::abs(norm_2(mLiftDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mLiftDirection);
            if (magnitude == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "lift_direction is not properly defined.",
                                   "")

                std::cout
                        << "WARNING: non unit vector detected in \"lift_direction\": "
                        << rParameters.PrettyPrintJsonString() << std::endl;
            std::cout << "normalizing \"lift_direction\"..." << std::endl;

            for (unsigned int d = 0; d < TDim; d++)
                mLiftDirection[d] /= magnitude;
        }

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~PotentialFlowLiftObjectiveFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if (rModelPart.HasSubModelPart(mStructureModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"structure_model_part_name\": ",
                mStructureModelPartName)
        }

// initialize the variables to zero.
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
                it->Set(STRUCTURE, false);
        }

        // mark structure
        ModelPart& rStructureModelPart = rModelPart.GetSubModelPart(mStructureModelPartName);
        for (auto it = rStructureModelPart.NodesBegin();
                it != rStructureModelPart.NodesEnd();
                ++it)
            it->Set(STRUCTURE, true);

        KRATOS_CATCH("")
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // allocate auxiliary memory. this is done here instead of Initialize()
        // in case of restart.
        int NumThreads = OpenMPUtils::GetNumThreads();
        mElementIds.resize(NumThreads);
        mDragFlagVector.resize(NumThreads);

        // use first element to initialize drag flag vector
        Element& rElem = *std::begin(rModelPart.Elements());
        #pragma omp parallel
        {
            // initialize drag flag and element id vectors
            int k = OpenMPUtils::ThisThread();
            mElementIds[k] = rElem.Id() + 1; // force initialization
            this->GetDragFlagVector(rElem);
        }

        KRATOS_CATCH("")
    }

    //here i do dL/dphi
    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
            const Matrix& rAdjointMatrix,
            Vector& rRHSContribution,
            ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        Vector& rDragFlagVector = this->GetDragFlagVector(rElem);
        noalias(rRHSContribution) = prod(rAdjointMatrix, rDragFlagVector);

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
            const Matrix& rAdjointMatrix,
            Vector& rRHSContribution,
            ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);
        rRHSContribution.clear();

        KRATOS_CATCH("")
    }

    
    //here i compute dL/dx
    virtual void CalculateSensitivityContribution(const Element& rElem,
            const Matrix& rDerivativesMatrix,
            Vector& rRHSContribution,
            ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);
        rRHSContribution.clear();


        Vector& rDragFlagVector = this->GetDragFlagVector(rElem);

        std::vector<unsigned int> active_nodes_index;
        active_nodes_index.reserve(2);

        const auto& geom = rElem.GetGeometry();
        for(unsigned int i=0; i<geom.size(); ++i)
            if(geom[i].Is(STRUCTURE))
                active_nodes_index.push_back(i);
            
        

        const unsigned int dim = rProcessInfo[DOMAIN_SIZE];
        if(dim == 2)
        {
            if(active_nodes_index.size() == 3)
                KRATOS_ERROR << "element " << rElem.Id() << " has 3 nodes on the boundary. can not compute contribution to Lift" << std::endl;
            else if(active_nodes_index.size() == 2) //do something ... otherwise return a vector of zeros
            {
                const array_1d<double,3> vinfinity = rProcessInfo[VELOCITY];
                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);
                
                bounded_matrix<double,3,2> DN_DX;
                array_1d<double,3> N;
                double vol;
                
                //calculate shape functions
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, vol);

                //gather nodal data
                Vector p(geom.size());
                for(unsigned int i=0; i<geom.size(); i++)
                {
                p[i] = geom[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                }

                const array_1d<double,2> v = prod(trans(DN_DX), p);

                const double rho = 1.0;
                const double pressure = 0.5*rho* (vinfinity_norm2 - inner_prod(v,v)); //0.5*(norm_2(vinfinity) - norm_2(v));
        
                //HERE WE COMPUTE pressure*diff(An,x);
                unsigned int i = active_nodes_index[0];
                unsigned int j = active_nodes_index[1];
                array_1d<double,2> x1,x2;
                x1[0] = geom[i].X();
                x1[1] = geom[i].Y();
                x2[0] = geom[j].X();
                x2[1] = geom[j].Y();

                array_1d<double,2> direction;
                direction[0] = mLiftDirection[0];
                direction[1] = mLiftDirection[1];
                
                //compute the pressure
                

                bounded_matrix<double,2,6> aux;
                aux.clear();
                aux(0,i*2+1) = 1.0*pressure;
                aux(0,j*2+1) = -1.0*pressure;
                aux(0,i*2) = -1.0*pressure;
                aux(0,j*2) = 1.0*pressure;

                noalias(rRHSContribution) = prod(direction,aux);



                //here we compute diff(pressure,x)*An
                array_1d<double,2> An;
                An[0] = -(x2[1]-x1[1]);
                An[1] = (x2[0]-x1[0]);

                bounded_matrix<double,TDim, TDim+1> x(TDim, TDim+1);
                for(unsigned int i=0; i<TDim+1; ++i)
                    for(unsigned int k=0; k<TDim; k++)
                        x(k,i) = this->GetGeometry()[i].Coordinates()[k];

                const double cdpAn_dx_out0 =             2*An[0];
                const double cdpAn_dx_out1 =             -x(2,0);
                const double cdpAn_dx_out2 =             cdpAn_dx_out1 + x(1,0);
                const double cdpAn_dx_out3 =             x(0,0) - x(1,0);
                const double cdpAn_dx_out4 =             cdpAn_dx_out3*p[2];
                const double cdpAn_dx_out5 =             cdpAn_dx_out1 + x(0,0);
                const double cdpAn_dx_out6 =             cdpAn_dx_out5*p[1];
                const double cdpAn_dx_out7 =             cdpAn_dx_out2*p[0] + cdpAn_dx_out4 - cdpAn_dx_out6;
                const double cdpAn_dx_out8 =             -x(2,1);
                const double cdpAn_dx_out9 =             cdpAn_dx_out8 + x(0,1);
                const double cdpAn_dx_out10 =             cdpAn_dx_out3*cdpAn_dx_out9;
                const double cdpAn_dx_out11 =             x(0,1) - x(1,1);
                const double cdpAn_dx_out12 =             cdpAn_dx_out11*cdpAn_dx_out5;
                const double cdpAn_dx_out13 =             cdpAn_dx_out10 - cdpAn_dx_out12;
                const double cdpAn_dx_out14 =             pow(cdpAn_dx_out13, -2);
                const double cdpAn_dx_out15 =             cdpAn_dx_out14*cdpAn_dx_out7;
                const double cdpAn_dx_out16 =             -p[2];
                const double cdpAn_dx_out17 =             cdpAn_dx_out8 + x(1,1);
                const double cdpAn_dx_out18 =             cdpAn_dx_out17*cdpAn_dx_out2*p[0];
                const double cdpAn_dx_out19 =             1.0/cdpAn_dx_out13;
                const double cdpAn_dx_out20 =             cdpAn_dx_out19*cdpAn_dx_out3*p[2];
                const double cdpAn_dx_out21 =             cdpAn_dx_out19*cdpAn_dx_out5*p[1];
                const double cdpAn_dx_out22 =             x(0,0)*x(1,1) - x(0,0)*x(2,1) - x(0,1)*x(1,0) + x(0,1)*x(2,0) + x(1,0)*x(2,1) - x(1,1)*x(2,0);
                const double cdpAn_dx_out23 =             pow(cdpAn_dx_out22, -2);
                const double cdpAn_dx_out24 =             cdpAn_dx_out17*cdpAn_dx_out23*p[0];
                const double cdpAn_dx_out25 =             cdpAn_dx_out11*cdpAn_dx_out14;
                const double cdpAn_dx_out26 =             cdpAn_dx_out25*p[2];
                const double cdpAn_dx_out27 =             cdpAn_dx_out14*p[1];
                const double cdpAn_dx_out28 =             cdpAn_dx_out27*cdpAn_dx_out9;
                const double cdpAn_dx_out29 =             p[0]/cdpAn_dx_out22;
                const double cdpAn_dx_out30 =             cdpAn_dx_out19*p[2];
                const double cdpAn_dx_out31 =             cdpAn_dx_out19*p[1];
                const double cdpAn_dx_out32 =             cdpAn_dx_out11*cdpAn_dx_out30 + cdpAn_dx_out17*cdpAn_dx_out29 - cdpAn_dx_out31*cdpAn_dx_out9;
                const double cdpAn_dx_out33 =             cdpAn_dx_out32*(cdpAn_dx_out24 + cdpAn_dx_out26 - cdpAn_dx_out28);
                const double cdpAn_dx_out34 =             cdpAn_dx_out15*(cdpAn_dx_out16 + cdpAn_dx_out17*cdpAn_dx_out20 - cdpAn_dx_out17*cdpAn_dx_out21 + cdpAn_dx_out18*cdpAn_dx_out19 + p[1]) + cdpAn_dx_out17*cdpAn_dx_out33;
                const double cdpAn_dx_out35 =             pow(cdpAn_dx_out7, 2)/pow(cdpAn_dx_out13, 3);
                const double cdpAn_dx_out36 =             cdpAn_dx_out2*cdpAn_dx_out35 + cdpAn_dx_out32*(cdpAn_dx_out18*cdpAn_dx_out23 + cdpAn_dx_out2*cdpAn_dx_out26 - cdpAn_dx_out2*cdpAn_dx_out28 + cdpAn_dx_out30 - cdpAn_dx_out31);
                const double cdpAn_dx_out37 =             cdpAn_dx_out19*cdpAn_dx_out2*p[0];
                const double cdpAn_dx_out38 =             cdpAn_dx_out15*(cdpAn_dx_out16 + cdpAn_dx_out20*cdpAn_dx_out9 - cdpAn_dx_out21*cdpAn_dx_out9 + cdpAn_dx_out37*cdpAn_dx_out9 + p[0]) + cdpAn_dx_out33*cdpAn_dx_out9;
                const double cdpAn_dx_out39 =             -cdpAn_dx_out29;
                const double cdpAn_dx_out40 =             cdpAn_dx_out32*(cdpAn_dx_out12*cdpAn_dx_out14*p[2] - cdpAn_dx_out14*cdpAn_dx_out6*cdpAn_dx_out9 + cdpAn_dx_out24*cdpAn_dx_out5 + cdpAn_dx_out30 + cdpAn_dx_out39) + cdpAn_dx_out35*cdpAn_dx_out5;
                const double cdpAn_dx_out41 =             cdpAn_dx_out11*cdpAn_dx_out33 + cdpAn_dx_out15*(cdpAn_dx_out11*cdpAn_dx_out20 - cdpAn_dx_out11*cdpAn_dx_out21 + cdpAn_dx_out11*cdpAn_dx_out37 + p[0] - p[1]);
                const double cdpAn_dx_out42 =             cdpAn_dx_out3*cdpAn_dx_out35 + cdpAn_dx_out32*(-cdpAn_dx_out10*cdpAn_dx_out27 + cdpAn_dx_out24*cdpAn_dx_out3 + cdpAn_dx_out25*cdpAn_dx_out4 + cdpAn_dx_out31 + cdpAn_dx_out39);
                const double cdpAn_dx_out43 =             2*An[1];
                aux(0,0)=-cdpAn_dx_out0*cdpAn_dx_out34;
                aux(0,1)=cdpAn_dx_out0*cdpAn_dx_out36;
                aux(0,2)=cdpAn_dx_out0*cdpAn_dx_out38;
                aux(0,3)=-cdpAn_dx_out0*cdpAn_dx_out40;
                aux(0,4)=-cdpAn_dx_out0*cdpAn_dx_out41;
                aux(0,5)=cdpAn_dx_out0*cdpAn_dx_out42;
                aux(1,0)=-cdpAn_dx_out34*cdpAn_dx_out43;
                aux(1,1)=cdpAn_dx_out36*cdpAn_dx_out43;
                aux(1,2)=cdpAn_dx_out38*cdpAn_dx_out43;
                aux(1,3)=-cdpAn_dx_out40*cdpAn_dx_out43;
                aux(1,4)=-cdpAn_dx_out41*cdpAn_dx_out43;
                aux(1,5)=cdpAn_dx_out42*cdpAn_dx_out43;

                noalias(rRHSContribution) += prod(direction,aux);


            }
            else
            {
                KRATOS_ERROR << " 3D CalculateSensitivityContribution is not yet implemented" << std::endl;
            }
        }

        KRATOS_CATCH("")
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mStructureModelPartName;
    array_1d<double, TDim> mLiftDirection;
    std::vector<Vector> mDragFlagVector;
    std::vector<unsigned int> mElementIds;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    Vector& GetDragFlagVector(const Element& rElement)
    {
        int k = OpenMPUtils::ThisThread();

        // if needed, compute the drag flag vector for this element
        if (rElement.Id() != mElementIds[k])
        {
            const unsigned int NumNodes = rElement.GetGeometry().PointsNumber();
            const unsigned int LocalSize = (TDim + 1) * NumNodes;

            if (mDragFlagVector[k].size() != LocalSize)
                mDragFlagVector[k].resize(LocalSize, false);

            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElement.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][LocalIndex++] = mLiftDirection[d];
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][LocalIndex++] = 0.0;
                }

                mDragFlagVector[k][LocalIndex++] = 0.0; // pressure dof
            }

            mElementIds[k] = rElement.Id();
        }

        return mDragFlagVector[k];
    }

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_POTENTIALFLOW_LIFT_OBJECTIVE_FUNCTION defined */
