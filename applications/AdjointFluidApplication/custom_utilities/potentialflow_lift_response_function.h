//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Iñigo Lopez and Riccardo Rossi
//

#if !defined(KRATOS_POTENTIALFLOW_LIFT_RESPONSE_FUNCTION)
#define KRATOS_POTENTIALFLOW_LIFT_RESPONSE_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/response_function.h"
#include "utilities/geometry_utilities.h"

//Other application includes
//#include "custom_conditions/potential_wall_condition.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for drag.
template <unsigned int TDim>
class PotentialFlowLiftResponseFunction : public ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PotentialFlowLiftResponseFunction);

    typedef ResponseFunction BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PotentialFlowLiftResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
    : ResponseFunction(rModelPart, rParameters)
    {
        KRATOS_TRY

        Parameters default_settings(R"(
        {
            "structure_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "lift_direction": [0.0, 1.0, 0.0]
        })");

        Parameters custom_settings = rParameters["custom_settings"];
        custom_settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName = custom_settings["structure_model_part_name"].GetString();
        

        if (custom_settings["lift_direction"].IsArray() == false ||
        custom_settings["lift_direction"].size() != 3)
        {
            KRATOS_ERROR << "lift_direction vector is not a vector or does "
                         <<  "not have size 3:"
                         <<  custom_settings.PrettyPrintJsonString();
        }

        for (unsigned int d = 0; d < 3; ++d)
            mLiftDirection[d] = custom_settings["lift_direction"][d].GetDouble();

        if (std::abs(norm_2(mLiftDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mLiftDirection);
            if (magnitude == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "lift_direction is not properly defined.",
                                   "")

                std::cout
                        << "WARNING: non unit vector detected in \"lift_direction\": "
                        << custom_settings.PrettyPrintJsonString() << std::endl;
            std::cout << "normalizing \"lift_direction\"..." << std::endl;

            for (unsigned int d = 0; d < 3; d++)
                mLiftDirection[d] /= magnitude;
        }

        //---------Added here by Inigo cause Initialized is not called------------------------------
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
            {
                it->Set(STRUCTURE, true);
            }
        //---------------------------------------------------------------------------------------

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~PotentialFlowLiftResponseFunction()
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
        std::cout << "initializing Response Function" << std::endl;

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
            {
                it->Set(STRUCTURE, true);
                it->Info();
            }
            

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

    /// Calculate the scalar valued response function
    double CalculateValue(ModelPart& rModelPart) override
    {
        //define and initialize lift variable
        double lift = 0;

        //getting processIngo
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        const array_1d<double,3> vinfinity = CurrentProcessInfo[VELOCITY];
        const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

        

        //loop over elements
        for (auto ielem = rModelPart.ElementsBegin();
                  ielem!= rModelPart.ElementsEnd();
                  ielem++)
        {
            const auto& geom = ielem->GetGeometry();    
            
            std::vector<unsigned int> active_nodes_index;
            active_nodes_index.reserve(2);//is this ok defined here?

            unsigned int counter = 0;
            for(unsigned int i=0; i<geom.size(); ++i)
            {
                if(geom[i].Is(STRUCTURE))
                {
                    active_nodes_index.push_back(i);
                    counter +=1;
                }
            }

            if(counter==TDim)//body element
            {
                if(ielem->IsNot(MARKER)) //normal element
                {
                    //compute Area*normal for 2D                     
                    unsigned int i = active_nodes_index[0];
                    unsigned int j = active_nodes_index[1];
                    
                    array_1d<double,2> x1,x2;
                    x1[0] = geom[i].X();
                    x1[1] = geom[i].Y();
                    x2[0] = geom[j].X();
                    x2[1] = geom[j].Y();
                    // std::cout << "X1 " << geom[0].X() << std::endl;
                    // std::cout << "Y1 " << geom[0].Y() << std::endl;
                    // std::cout << "X2 " << geom[1].X() << std::endl;
                    // std::cout << "Y2 " << geom[1].Y() << std::endl;
                    // std::cout << "X3 " << geom[2].X() << std::endl;
                    // std::cout << "Y3 " << geom[2].Y() << std::endl;
                              
                    array_1d<double,3> An;
                    An[0] = (x2[1]-x1[1]);//carefull with signs! normal pointing outwards in this case
                    An[1] = -(x2[0]-x1[0]);
                    An[2] = 0; //this is only valid for 2D!!
                    
                    //compute pressure coefficient cp
                    bounded_matrix<double,3,2> DN_DX;
                    array_1d<double,3> N;
                    double vol;
                    
                    //calculate shape functions
                    GeometryUtils::CalculateGeometryData(ielem->GetGeometry(), DN_DX, N, vol);
                    //gather nodal data
                    Vector p(geom.size());
                    for(unsigned int i=0; i<geom.size(); i++)
                    {
                        p[i] = geom[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                    }
                    const array_1d<double,2> v = prod(trans(DN_DX), p);
                    //const double rho = 1.0;
                    const double cp = (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*rho* (vinfinity_norm2 - inner_prod(v,v)); //
                         //compute lift
                    lift += cp*inner_prod(An,mLiftDirection); 
    
                }
                else //wake-body element (only one in 2d)
                {
                    //compute Area*normal for 2D
                    unsigned int i = active_nodes_index[0];
                    unsigned int j = active_nodes_index[1];
                    
                    array_1d<double,2> x1,x2;
                    x1[0] = geom[i].X();
                    x1[1] = geom[i].Y();
                    x2[0] = geom[j].X();
                    x2[1] = geom[j].Y();
                              
                    array_1d<double,3> An;
                    An[0] = (x2[1]-x1[1]);//carefull with signs! normal pointing outwards in this case
                    An[1] = -(x2[0]-x1[0]);
                    An[2] = 0; //this is only valid for 2D!!
                    
                    //compute pressure coefficient cp
                    bounded_matrix<double,3,2> DN_DX;
                    array_1d<double,3> N;
                    double vol;
                    
                    //calculate shape functions
                    GeometryUtils::CalculateGeometryData(ielem->GetGeometry(), DN_DX, N, vol);
                    
                    //gather nodal data
                    Vector p(geom.size());
                    Vector Potential = ZeroVector(2*geom.size());
                    ielem->GetFirstDerivativesVector(Potential,0);//this takes both upper and lower potential
                    for(unsigned int i=0; i<geom.size(); i++)//extracting only the positive side
                    {
                        p[i] = Potential[i];
                    }

                    const array_1d<double,2> v = prod(trans(DN_DX), p);
                    //const double rho = 1.0;
                    const double cp = (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*rho* (vinfinity_norm2 - inner_prod(v,v)); //
                    //compute lift
                    lift += cp*inner_prod(An,mLiftDirection);    
                }

            }
            else if(counter>TDim)
            {
                KRATOS_ERROR << "element " << ielem->Id() << " has " << counter << "nodes on the boundary. cannot compute contribution to Lift" << std::endl;                
            }
            else//counter = 1 or counter = 0 not an airfoil element
            {
                lift +=0;
            }                       
        }

        // //loop over conditions
        // for (auto icond = rModelPart.ConditionsBegin();
        //          icond!=rModelPart.ConditionsEnd(); icond++)
        // {
        //     //compute Area*normal
        //     array_1d<double,3> An;

        //     icond->Calculate(ADJOINT_VELOCITY, An, CurrentProcessInfo);//TODO: there has to be a better way of doing this

        //     std::cout << "An " << An << std::endl;

        //     //get pressure
        //     double pressure = icond->GetValue(PRESSURE);

        //     //compute lift
        //     lift += pressure*inner_prod(An,mLiftDirection);
        // }
        return lift;
    }

    /// Calculate the lift local gradient w.r.t. primal solution (i.e. potential)
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. primal.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
     void CalculateGradient(Element& rAdjointElem,
        const Matrix& rAdjointMatrix,
        Vector& rResponseGradient,
        ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        const auto& geom = rAdjointElem.GetGeometry();

        std::vector<unsigned int> active_nodes_index;
        active_nodes_index.reserve(2);        
                
        unsigned int counter = 0;
        for(unsigned int i=0; i<geom.size(); ++i)
        if(geom[i].Is(STRUCTURE))
        {
            active_nodes_index.push_back(i);
            counter +=1;
        }

        if(counter==TDim)//body element
        {
            if(rAdjointElem.IsNot(MARKER)) //normal element
            {
                //std::cout << "NORMAL ELEMENT ENTERING CalculateGradient from potentialflow_lift_response_function" << std::endl;
                //making sure rResponseGradient has the correct size
                if (rResponseGradient.size() != rAdjointMatrix.size1())
                rResponseGradient.resize(rAdjointMatrix.size1(), false);
                rResponseGradient.clear();

                const array_1d<double,3> vinfinity = rProcessInfo[VELOCITY];
                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                //dividing the LeftHandSideMatrix by the element's volume/area
                Matrix LeftHandSideMatrix = trans(rAdjointMatrix);
                if(TDim == 2) LeftHandSideMatrix /= rAdjointElem.GetGeometry().Area();
                else LeftHandSideMatrix /= rAdjointElem.GetGeometry().Volume();

                //computing the residual (wo volume/area)
                Vector Potential = ZeroVector(rAdjointMatrix.size1());
                rAdjointElem.GetFirstDerivativesVector(Potential,0);
                auto Residual = -prod(LeftHandSideMatrix,Potential)/vinfinity_norm2;//TODO: check sign
                //std::cout << "Residual = " << Residual << std::endl;

                //getting Area*normal from the condition
                // array_1d<double,3> An;
                // rAdjointElem.Calculate(ADJOINT_VELOCITY,An,rProcessInfo);//this needs to be implemented
                //compute Area*normal for 2D
                
                unsigned int i = active_nodes_index[0];
                unsigned int j = active_nodes_index[1];

                //std::cout << "active_nodes_index[0] =" << active_nodes_index[0] << std::endl;
                //std::cout << "active_nodes_index[1] =" << active_nodes_index[1] << std::endl;
                array_1d<double,2> x1,x2;
                x1[0] = geom[i].X();
                x1[1] = geom[i].Y();
                x2[0] = geom[j].X();
                x2[1] = geom[j].Y();

                array_1d<double,3> An;
                An[0] = (x2[1]-x1[1]);//carefull with signs! normal pointing outwards in this case
                An[1] = -(x2[0]-x1[0]);
                An[2] = 0; //this is only valid for 2D!!

                //compute dL/dPhi
                rResponseGradient = 2*Residual*inner_prod(An,mLiftDirection);
                //std::cout << "rResponseGradient =" << rResponseGradient << std::endl;
                
                //std::cout << "NORMAL ELEMENT EXITING CalculateGradient from potentialflow_lift_response_function" << std::endl;
            }
            else //wake-body element (only one in 2d)
            {
                //std::cout << "WAKE-BODY ELEMENT ENTERING CalculateGradient from potentialflow_lift_response_function" << std::endl;
                //making sure rResponseGradient has the correct size (6)
                if (rResponseGradient.size() != rAdjointMatrix.size1())
                rResponseGradient.resize(rAdjointMatrix.size1(), false);
                rResponseGradient.clear();

                const array_1d<double,3> vinfinity = rProcessInfo[VELOCITY];
                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                bounded_matrix<double,3,2> DN;
                array_1d<double,3> N;
                double vol;
                
                //calculate shape functions
                GeometryUtils::CalculateGeometryData(rAdjointElem.GetGeometry(), DN, N, vol);

                //gather nodal data
                Vector p(geom.size());
                Vector Potential = ZeroVector(2*geom.size());
                rAdjointElem.GetFirstDerivativesVector(Potential,0);//this takes both upper and lower potential
                for(unsigned int i=0; i<geom.size(); i++)//extracting only the negative (lower) side
                {
                    p[i] = Potential[geom.size()+i];
                }

                //computing dcp/dPhi, which is equivalent to the residual (wo volume/area)
                //auto dcpdPhi = -prod(prod(DN,trans(DN)),p)/vinfinity_norm2;

                //getting Area*normal from the condition
                // array_1d<double,3> An;
                // rAdjointElem.Calculate(ADJOINT_VELOCITY,An,rProcessInfo);//this needs to be implemented
                //compute Area*normal for 2D
                
                unsigned int i = active_nodes_index[0];
                unsigned int j = active_nodes_index[1];

                //std::cout << "active_nodes_index[0] =" << active_nodes_index[0] << std::endl;
                //std::cout << "active_nodes_index[1] =" << active_nodes_index[1] << std::endl;
                array_1d<double,2> x1,x2;
                x1[0] = geom[i].X();
                x1[1] = geom[i].Y();
                x2[0] = geom[j].X();
                x2[1] = geom[j].Y();

                array_1d<double,3> An;
                An[0] = (x2[1]-x1[1]);//carefull with signs! normal pointing outwards in this case
                An[1] = -(x2[0]-x1[0]);
                An[2] = 0; //this is only valid for 2D!!

                auto dcldPhi = -2*prod(prod(DN,trans(DN)),p)*inner_prod(An,mLiftDirection)/vinfinity_norm2;

                //compute dL/dPhi
                for(unsigned int i=0; i<geom.size(); i++)//extracting only the negative (lower) side
                {
                    rResponseGradient[i] = 0;// dcldPhi(i);
                    //rResponseGradient[i+geom.size()] = 2*dcpdPhi[i]*inner_prod(An,mLiftDirection);
                    rResponseGradient[i+geom.size()] = dcldPhi(i);
                }
                //rResponseGradient = 2*dcpdPhi*inner_prod(An,mLiftDirection);
                //std::cout << "rResponseGradient =" << rResponseGradient << std::endl;
                
                //std::cout << "WAKE-BODY ELEMENT EXITING CalculateGradient from potentialflow_lift_response_function" << std::endl;
            }
        }
        else if(counter>TDim)
        {
            KRATOS_ERROR << "element " << rAdjointElem.Id() << " has " << counter << "nodes on the boundary. Cannot compute contribution to Lift" << std::endl;                
        }
        else//counter = 1 or counter = 0 not an airfoil element
        {
            //making sure rResponseGradient has the correct size
            if (rResponseGradient.size() != rAdjointMatrix.size1())
                rResponseGradient.resize(rAdjointMatrix.size1(), false);
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    //here i do dL/dphi
    // virtual void CalculateAdjointVelocityContribution(const Element& rElem,
    //         const Matrix& rAdjointMatrix,
    //         Vector& rRHSContribution,
    //         ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY

    //     if (rRHSContribution.size() != rAdjointMatrix.size1())
    //         rRHSContribution.resize(rAdjointMatrix.size1(), false);

    //     Vector& rDragFlagVector = this->GetDragFlagVector(rElem);
    //     noalias(rRHSContribution) = prod(rAdjointMatrix, rDragFlagVector);

    //     KRATOS_CATCH("")
    // }

    // virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
    //         const Matrix& rAdjointMatrix,
    //         Vector& rRHSContribution,
    //         ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY

    //     if (rRHSContribution.size() != rAdjointMatrix.size1())
    //         rRHSContribution.resize(rAdjointMatrix.size1(), false);
    //     rRHSContribution.clear();

    //     KRATOS_CATCH("")
    // }

    
    //here i compute dL/dx
    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
     void CalculateSensitivityGradient(Element& rAdjointElem,
                                       const Variable<array_1d<double,3>>& rVariable,
                                       const Matrix& rDerivativesMatrix,
                                       Vector& rResponseGradient,
                                       ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
        rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();


        //Vector& rDragFlagVector = this->GetDragFlagVector(rAdjointElem);

        std::vector<unsigned int> active_nodes_index;
        active_nodes_index.reserve(2);

        const auto& geom = rAdjointElem.GetGeometry();
        int counter = 0;
        for(unsigned int i=0; i<geom.size(); ++i)
            if(geom[i].Is(STRUCTURE))
            {
                active_nodes_index.push_back(i);
                counter +=1;                
            }

        //const unsigned int dim = rProcessInfo[DOMAIN_SIZE];
        if(TDim == 2)
        {
            if(counter == 3)
                KRATOS_ERROR << "element " << rAdjointElem.Id() << " has 3 nodes on the boundary. cannot compute contribution to Lift" << std::endl;            
            else if(counter == 2) //do something ... otherwise return a vector of zeros
            {
                if(rAdjointElem.IsNot(MARKER))//normal element
                {
                    // int kutta_check = 0;
                    // for(unsigned int i=0; i<geom.size(); i++)
                    // {
                    //     if(fabs(geom[i].X() - 1.0) < 0.02)
                    //     {
                    //         //std::cout << "I AM A KUTTA ELEMENT!" << std::endl;
                    //         kutta_check = 1;
                    //     } 
                    // }

                    const array_1d<double,3> vinfinity = rProcessInfo[VELOCITY];
                    const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);
                    
                    bounded_matrix<double,3,2> DN;
                    array_1d<double,3> N;
                    double vol;
                    
                    //calculate shape functions
                    GeometryUtils::CalculateGeometryData(rAdjointElem.GetGeometry(), DN, N, vol);
    
                    //gather nodal data
                    Vector p(geom.size());
                    for(unsigned int i=0; i<geom.size(); i++)
                    {
                        p[i] = geom[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                    }
    
                    const array_1d<double,2> v = prod(trans(DN), p);
    
                    //const double rho = 1.0;
                    //compute the pressure
                    const double pressure = (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*rho* (vinfinity_norm2 - inner_prod(v,v)); //
            
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
                    
                    bounded_matrix<double,2,6> aux;
                    aux.clear();
                    aux(0,i*2+1) = -1.0*pressure;
                    aux(0,j*2+1) = 1.0*pressure;
                    aux(1,i*2) = 1.0*pressure;
                    aux(1,j*2) = -1.0*pressure;
    
                    noalias(rResponseGradient) = prod(direction,aux);
                    //std::cout << "rResponseGradient =" << rResponseGradient << std::endl;
    
                    //here we compute diff(pressure,x)*An
                    array_1d<double,2> An;
                    An[0] = (x2[1]-x1[1]);//carfull with signs! normal pointing outwards in this case
                    An[1] = -(x2[0]-x1[0]);
    
                    //const unsigned int nnodes = rAdjointElem.GetGeometry().size();
    
                    // //matrix of coordinates
                    // bounded_matrix<double,TDim, nnodes> x(nnodes,TDim);
                    // for(unsigned int i=0; i<nnodes; ++i)
                    //     for(unsigned int k=0; k<TDim; k++)
                    //         x(i,k) = rAdjointElem.GetGeometry()[i].Coordinates()[k];
    
                    bounded_matrix<double,TDim, TDim+1> x(TDim+1, TDim);
                    for(unsigned int i=0; i<TDim+1; ++i)
                        for(unsigned int k=0; k<TDim; k++)
                            x(i,k) = rAdjointElem.GetGeometry()[i].Coordinates()[k];
    
                    // bounded_matrix<double,TDim, TDim+1> x(TDim, TDim+1);
                    // for(unsigned int i=0; i<TDim+1; ++i)
                    //     for(unsigned int k=0; k<TDim; k++)
                    //         x(k,i) = rAdjointElem.GetGeometry()[i].Coordinates()[k];
    
    
                    const double caux0 =             1.0/vinfinity_norm2;
                    const double caux1 =             x(0,0) - x(1,0);
                    const double caux2 =             -x(2,1);
                    const double caux3 =             caux2 + x(0,1);
                    const double caux4 =             -x(2,0);
                    const double caux5 =             caux4 + x(0,0);
                    const double caux6 =             x(0,1) - x(1,1);
                    const double caux7 =             caux1*caux3 - caux5*caux6;
                    const double caux8 =             pow(caux7, -3);
                    const double caux9 =             2*An[0]*caux0*caux8;
                    const double caux10 =             caux4 + x(1,0);
                    const double caux11 =             caux1*p[2] + caux10*p[0] - caux5*p[1];
                    const double caux12 =             -p[2];
                    const double caux13 =             caux7*(caux12 + p[1]);
                    const double caux14 =             caux2 + x(1,1);
                    const double caux15 =             caux14*p[0] - caux3*p[1] + caux6*p[2];
                    const double caux16 =             pow(caux11, 2) + pow(caux15, 2);
                    const double caux17 =             caux11*caux13 + caux14*caux16;
                    const double caux18 =             -caux10*caux16 + caux13*caux15;
                    const double caux19 =             caux7*(caux12 + p[0]);
                    const double caux20 =             caux11*caux19 + caux16*caux3;
                    const double caux21 =             caux15*caux19 - caux16*caux5;
                    const double caux22 =             caux7*(p[0] - p[1]);
                    const double caux23 =             caux11*caux22 + caux16*caux6;
                    const double caux24 =             -caux1*caux16 + caux15*caux22;
                    const double caux25 =             2*An[1]*caux0*caux8;
                    aux(0,0)=caux17*caux9;
                    aux(0,1)=caux18*caux9;
                    aux(0,2)=-caux20*caux9;
                    aux(0,3)=-caux21*caux9;
                    aux(0,4)=caux23*caux9;
                    aux(0,5)=caux24*caux9;
                    aux(1,0)=caux17*caux25;
                    aux(1,1)=caux18*caux25;
                    aux(1,2)=-caux20*caux25;
                    aux(1,3)=-caux21*caux25;
                    aux(1,4)=caux23*caux25;
                    aux(1,5)=caux24*caux25;                
    
                    // const double caux0 =             1.0/vinfinity_norm2;
                    // const double caux1 =             DN(0,0)*x(0,0) + DN(1,0)*x(1,0) + DN(2,0)*x(2,0);
                    // const double caux2 =             DN(0,1)*x(0,1) + DN(1,1)*x(1,1) + DN(2,1)*x(2,1);
                    // const double caux3 =             DN(0,0)*x(0,1) + DN(1,0)*x(1,1) + DN(2,0)*x(2,1);
                    // const double caux4 =             DN(0,1)*x(0,0) + DN(1,1)*x(1,0) + DN(2,1)*x(2,0);
                    // const double caux5 =             caux1*caux2 - caux3*caux4;
                    // const double caux6 =             pow(caux5, -3);
                    // const double caux7 =             2*An[0]*caux0*caux6;
                    // const double caux8 =             DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0);
                    // const double caux9 =             DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0);
                    // const double caux10 =             caux5*(caux8*p[1] + caux9*p[2]);
                    // const double caux11 =             DN(0,0)*caux4 - DN(0,1)*caux1;
                    // const double caux12 =             DN(1,0)*caux4 - DN(1,1)*caux1;
                    // const double caux13 =             DN(2,0)*caux4 - DN(2,1)*caux1;
                    // const double caux14 =             caux11*p[0] + caux12*p[1] + caux13*p[2];
                    // const double caux15 =             DN(0,0)*caux2 - DN(0,1)*caux3;
                    // const double caux16 =             DN(1,0)*caux2 - DN(1,1)*caux3;
                    // const double caux17 =             DN(2,0)*caux2 - DN(2,1)*caux3;
                    // const double caux18 =             caux15*p[0] + caux16*p[1] + caux17*p[2];
                    // const double caux19 =             pow(caux14, 2) + pow(caux18, 2);
                    // const double caux20 =             caux10*caux14 + caux15*caux19;
                    // const double caux21 =             -caux10*caux18 + caux11*caux19;
                    // const double caux22 =             DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0);
                    // const double caux23 =             caux5*(-caux22*p[2] + caux8*p[0]);
                    // const double caux24 =             -caux14*caux23 + caux16*caux19;
                    // const double caux25 =             caux12*caux19 + caux18*caux23;
                    // const double caux26 =             caux5*(caux22*p[1] + caux9*p[0]);
                    // const double caux27 =             -caux14*caux26 + caux17*caux19;
                    // const double caux28 =             caux13*caux19 + caux18*caux26;
                    // const double caux29 =             2*An[1]*caux0*caux6;
                    // aux(0,0)=caux20*caux7;
                    // aux(0,1)=-caux21*caux7;
                    // aux(0,2)=caux24*caux7;
                    // aux(0,3)=-caux25*caux7;
                    // aux(0,4)=caux27*caux7;
                    // aux(0,5)=-caux28*caux7;
                    // aux(1,0)=caux20*caux29;
                    // aux(1,1)=-caux21*caux29;
                    // aux(1,2)=caux24*caux29;
                    // aux(1,3)=-caux25*caux29;
                    // aux(1,4)=caux27*caux29;
                    // aux(1,5)=-caux28*caux29;
                                            
                    
                    // const double cdpAn_dx_out0 =             2*An[0];
                    // const double cdpAn_dx_out1 =             -x(2,0);
                    // const double cdpAn_dx_out2 =             cdpAn_dx_out1 + x(1,0);
                    // const double cdpAn_dx_out3 =             x(0,0) - x(1,0);
                    // const double cdpAn_dx_out4 =             cdpAn_dx_out3*p[2];
                    // const double cdpAn_dx_out5 =             cdpAn_dx_out1 + x(0,0);
                    // const double cdpAn_dx_out6 =             cdpAn_dx_out5*p[1];
                    // const double cdpAn_dx_out7 =             cdpAn_dx_out2*p[0] + cdpAn_dx_out4 - cdpAn_dx_out6;
                    // const double cdpAn_dx_out8 =             -x(2,1);
                    // const double cdpAn_dx_out9 =             cdpAn_dx_out8 + x(0,1);
                    // const double cdpAn_dx_out10 =             cdpAn_dx_out3*cdpAn_dx_out9;
                    // const double cdpAn_dx_out11 =             x(0,1) - x(1,1);
                    // const double cdpAn_dx_out12 =             cdpAn_dx_out11*cdpAn_dx_out5;
                    // const double cdpAn_dx_out13 =             cdpAn_dx_out10 - cdpAn_dx_out12;
                    // const double cdpAn_dx_out14 =             pow(cdpAn_dx_out13, -2);
                    // const double cdpAn_dx_out15 =             cdpAn_dx_out14*cdpAn_dx_out7;
                    // const double cdpAn_dx_out16 =             -p[2];
                    // const double cdpAn_dx_out17 =             cdpAn_dx_out8 + x(1,1);
                    // const double cdpAn_dx_out18 =             cdpAn_dx_out17*cdpAn_dx_out2*p[0];
                    // const double cdpAn_dx_out19 =             1.0/cdpAn_dx_out13;
                    // const double cdpAn_dx_out20 =             cdpAn_dx_out19*cdpAn_dx_out3*p[2];
                    // const double cdpAn_dx_out21 =             cdpAn_dx_out19*cdpAn_dx_out5*p[1];
                    // const double cdpAn_dx_out22 =             x(0,0)*x(1,1) - x(0,0)*x(2,1) - x(0,1)*x(1,0) + x(0,1)*x(2,0) + x(1,0)*x(2,1) - x(1,1)*x(2,0);
                    // const double cdpAn_dx_out23 =             pow(cdpAn_dx_out22, -2);
                    // const double cdpAn_dx_out24 =             cdpAn_dx_out17*cdpAn_dx_out23*p[0];
                    // const double cdpAn_dx_out25 =             cdpAn_dx_out11*cdpAn_dx_out14;
                    // const double cdpAn_dx_out26 =             cdpAn_dx_out25*p[2];
                    // const double cdpAn_dx_out27 =             cdpAn_dx_out14*p[1];
                    // const double cdpAn_dx_out28 =             cdpAn_dx_out27*cdpAn_dx_out9;
                    // const double cdpAn_dx_out29 =             p[0]/cdpAn_dx_out22;
                    // const double cdpAn_dx_out30 =             cdpAn_dx_out19*p[2];
                    // const double cdpAn_dx_out31 =             cdpAn_dx_out19*p[1];
                    // const double cdpAn_dx_out32 =             cdpAn_dx_out11*cdpAn_dx_out30 + cdpAn_dx_out17*cdpAn_dx_out29 - cdpAn_dx_out31*cdpAn_dx_out9;
                    // const double cdpAn_dx_out33 =             cdpAn_dx_out32*(cdpAn_dx_out24 + cdpAn_dx_out26 - cdpAn_dx_out28);
                    // const double cdpAn_dx_out34 =             cdpAn_dx_out15*(cdpAn_dx_out16 + cdpAn_dx_out17*cdpAn_dx_out20 - cdpAn_dx_out17*cdpAn_dx_out21 + cdpAn_dx_out18*cdpAn_dx_out19 + p[1]) + cdpAn_dx_out17*cdpAn_dx_out33;
                    // const double cdpAn_dx_out35 =             pow(cdpAn_dx_out7, 2)/pow(cdpAn_dx_out13, 3);
                    // const double cdpAn_dx_out36 =             cdpAn_dx_out2*cdpAn_dx_out35 + cdpAn_dx_out32*(cdpAn_dx_out18*cdpAn_dx_out23 + cdpAn_dx_out2*cdpAn_dx_out26 - cdpAn_dx_out2*cdpAn_dx_out28 + cdpAn_dx_out30 - cdpAn_dx_out31);
                    // const double cdpAn_dx_out37 =             cdpAn_dx_out19*cdpAn_dx_out2*p[0];
                    // const double cdpAn_dx_out38 =             cdpAn_dx_out15*(cdpAn_dx_out16 + cdpAn_dx_out20*cdpAn_dx_out9 - cdpAn_dx_out21*cdpAn_dx_out9 + cdpAn_dx_out37*cdpAn_dx_out9 + p[0]) + cdpAn_dx_out33*cdpAn_dx_out9;
                    // const double cdpAn_dx_out39 =             -cdpAn_dx_out29;
                    // const double cdpAn_dx_out40 =             cdpAn_dx_out32*(cdpAn_dx_out12*cdpAn_dx_out14*p[2] - cdpAn_dx_out14*cdpAn_dx_out6*cdpAn_dx_out9 + cdpAn_dx_out24*cdpAn_dx_out5 + cdpAn_dx_out30 + cdpAn_dx_out39) + cdpAn_dx_out35*cdpAn_dx_out5;
                    // const double cdpAn_dx_out41 =             cdpAn_dx_out11*cdpAn_dx_out33 + cdpAn_dx_out15*(cdpAn_dx_out11*cdpAn_dx_out20 - cdpAn_dx_out11*cdpAn_dx_out21 + cdpAn_dx_out11*cdpAn_dx_out37 + p[0] - p[1]);
                    // const double cdpAn_dx_out42 =             cdpAn_dx_out3*cdpAn_dx_out35 + cdpAn_dx_out32*(-cdpAn_dx_out10*cdpAn_dx_out27 + cdpAn_dx_out24*cdpAn_dx_out3 + cdpAn_dx_out25*cdpAn_dx_out4 + cdpAn_dx_out31 + cdpAn_dx_out39);
                    // const double cdpAn_dx_out43 =             2*An[1];
                    // aux(0,0)=-cdpAn_dx_out0*cdpAn_dx_out34;
                    // aux(0,1)=cdpAn_dx_out0*cdpAn_dx_out36;
                    // aux(0,2)=cdpAn_dx_out0*cdpAn_dx_out38;
                    // aux(0,3)=-cdpAn_dx_out0*cdpAn_dx_out40;
                    // aux(0,4)=-cdpAn_dx_out0*cdpAn_dx_out41;
                    // aux(0,5)=cdpAn_dx_out0*cdpAn_dx_out42;
                    // aux(1,0)=-cdpAn_dx_out34*cdpAn_dx_out43;
                    // aux(1,1)=cdpAn_dx_out36*cdpAn_dx_out43;
                    // aux(1,2)=cdpAn_dx_out38*cdpAn_dx_out43;
                    // aux(1,3)=-cdpAn_dx_out40*cdpAn_dx_out43;
                    // aux(1,4)=-cdpAn_dx_out41*cdpAn_dx_out43;
                    // aux(1,5)=cdpAn_dx_out42*cdpAn_dx_out43;
    
                    // std::cout << "caux0   =" << caux0  << std::endl;            
                    // std::cout << "caux1   =" << caux1  << std::endl;
                    // std::cout << "caux2   =" << caux2  << std::endl;
                    // std::cout << "caux3   =" << caux3  << std::endl;
                    // std::cout << "caux4   =" << caux4  << std::endl;
                    // std::cout << "caux5   =" << caux5  << std::endl;
                    // std::cout << "caux6   =" << caux6  << std::endl;
                    // std::cout << "caux7   =" << caux7  << std::endl;
                    // std::cout << "caux8   =" << caux8  << std::endl;
                    // std::cout << "caux9   =" << caux9  << std::endl;
                    // std::cout << "caux10  =" << caux10 << std::endl;   
                    // std::cout << "caux11  =" << caux11 << std::endl;
                    // std::cout << "caux12  =" << caux12 << std::endl;
                    // std::cout << "caux13  =" << caux13 << std::endl;
                    // std::cout << "caux14  =" << caux14 << std::endl;
                    // std::cout << "caux15  =" << caux15 << std::endl;
                    // std::cout << "caux16  =" << caux16 << std::endl;
                    // std::cout << "caux17  =" << caux17 << std::endl;
                    // std::cout << "caux18  =" << caux18 << std::endl;
                    // std::cout << "caux19  =" << caux19 << std::endl;
                    // std::cout << "caux20  =" << caux20 << std::endl;
                    // std::cout << "caux21  =" << caux21 << std::endl;
                    // std::cout << "caux22  =" << caux22 << std::endl;
                    // std::cout << "caux23  =" << caux23 << std::endl;
                    // std::cout << "caux24  =" << caux24 << std::endl;
                    // std::cout << "caux25  =" << caux25 << std::endl;
                    // std::cout << "caux26  =" << caux26 << std::endl;
                    // std::cout << "caux27  =" << caux27 << std::endl;
                    // std::cout << "caux28  =" << caux28 << std::endl;
                    // std::cout << "caux29  =" << caux29 << std::endl;
                    // std::cout << "cdpAn_dx_out30    =" << cdpAn_dx_out30 << std::endl;
                    // std::cout << "cdpAn_dx_out31    =" << cdpAn_dx_out31 << std::endl;
                    // std::cout << "cdpAn_dx_out32    =" << cdpAn_dx_out32 << std::endl;
                    // std::cout << "cdpAn_dx_out33    =" << cdpAn_dx_out33 << std::endl;
                    // std::cout << "cdpAn_dx_out34    =" << cdpAn_dx_out34 << std::endl;
                    // std::cout << "cdpAn_dx_out35    =" << cdpAn_dx_out35 << std::endl;
                    // std::cout << "cdpAn_dx_out36    =" << cdpAn_dx_out36 << std::endl;
                    // std::cout << "cdpAn_dx_out37    =" << cdpAn_dx_out37 << std::endl;
                    // std::cout << "cdpAn_dx_out38    =" << cdpAn_dx_out38 << std::endl;
                    // std::cout << "cdpAn_dx_out39    =" << cdpAn_dx_out39 << std::endl;
                    // std::cout << "cdpAn_dx_out40    =" << cdpAn_dx_out40 << std::endl;
                    // std::cout << "cdpAn_dx_out41    =" << cdpAn_dx_out41 << std::endl;
                    // std::cout << "cdpAn_dx_out42    =" << cdpAn_dx_out42 << std::endl;
                    // std::cout << "cdpAn_dx_out43    =" << cdpAn_dx_out43 << std::endl;             
                    
                    //std::cout << "aux =" << aux << std::endl;
    
                    noalias(rResponseGradient) += prod(direction,aux);
                    //std::cout << "SensitivityGradient dL/dx =" << rResponseGradient << std::endl;
                    //const double sensitivity_norm2 = inner_prod(rResponseGradient,rResponseGradient);
                    // if(rResponseGradient(0)>100||rResponseGradient(1)>100||rResponseGradient(2)>100)
                    // {
                    //     std::cout << "I AM HERE! dL/dx ="<< std::endl;
                    //     std::cout << "SensitivityGradient dL/dx =" << rResponseGradient << std::endl;
                    //     std::cout << "sensitivity_norm2 dL/dx =" << sensitivity_norm2 << std::endl;
                    //     std::cout << "geom =" << geom << std::endl;
                    // }
                    // if(sensitivity_norm2>9)
                    // {
                    //     std::cout << "SensitivityGradient dL/dx =" << rResponseGradient << std::endl;
                    //     std::cout << "sensitivity_norm2 dL/dx =" << sensitivity_norm2 << std::endl;
                    //     std::cout << "geom =" << geom << std::endl;
                    // }
                    // if(kutta_check==1)
                    // {
                    //     std::cout << "I AM A KUTTA ELEMENT!" << std::endl;
                    //     rResponseGradient.clear(); 
                    // }

                }
                else//wake element
                {
                    std::cout << "I AM THE ONLY STRUCTURE-WAKE CASE" << std::endl;
                    rResponseGradient.clear();
                    // std::cout << "SensitivityGradient dL/dx =" << rResponseGradient << std::endl;
                    // std::cout << "geom =" << geom << std::endl;
                }
            }
            else//counter = 1 or counter = 0 not an airfoil element
            {
                // std::cout << "I AM ENTERING A Non-airfoil element" << std::endl;
                // std::cout << "counter =" << counter << std::endl;
                rResponseGradient.clear();
            }                
        }
        else if(TDim == 3)
        {
            KRATOS_ERROR << " 3D CalculateSensitivityContribution is not yet implemented" << std::endl;
        }
        else
        {
            KRATOS_ERROR << " Wrong Dimentsion" << std::endl;
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
    array_1d<double, 3> mLiftDirection;
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

#endif /* KRATOS_POTENTIALFLOW_LIFT_RESPONSE_FUNCTION defined */
