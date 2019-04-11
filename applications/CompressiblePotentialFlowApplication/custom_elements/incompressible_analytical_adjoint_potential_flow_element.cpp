//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Iñigo Lopez based on M. Nuñez, A. Geiser, M. Fusseder and R. Rossi work
//
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "incompressible_analytical_adjoint_potential_flow_element.h"
#include "custom_elements/potential_flow_functions.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
          return Kratos::make_shared<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
            return Kratos::make_shared<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const 
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointAnalyticalIncompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::IntegrationMethod AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetIntegrationMethod() const 
    {
        return mpPrimalElement->GetIntegrationMethod();
    }
  
    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Initialize() 
    {   
        mpPrimalElement->Initialize();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->Data() = this->Data();
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {

    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) 
    {
        MatrixType tmp;
        mpPrimalElement->CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);                                    
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
    
    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) 
    {
        KRATOS_TRY
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if (wake == 1) // wake element
        {
            if(rValues.size() != 2*NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);
            GetValuesOnSplitElement(rValues,distances);
            
        }else{ // normal element
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if(kutta == 0){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] =GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY;
        const double delta = this->GetPerturbationSize();
        ProcessInfo process_info = rCurrentProcessInfo;

        Vector RHS;
        Vector RHS_perturbed;

        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);

        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                if ((GetGeometry()[i_node].Is(SOLID)) && (!GetGeometry()[i_node].GetValue(TRAILING_EDGE))){
                    pGetPrimalElement()->GetGeometry()[i_node].GetInitialPosition()[i_dim] += delta;
                    pGetPrimalElement()->GetGeometry()[i_node].Coordinates()[i_dim] += delta;

                    // compute LHS after perturbation
                    pGetPrimalElement()->CalculateRightHandSide(RHS_perturbed, process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = (RHS_perturbed[i] - RHS[i]) / delta;

                    // unperturb the design variable
                    pGetPrimalElement()->GetGeometry()[i_node].GetInitialPosition()[i_dim] -= delta;
                    pGetPrimalElement()->GetGeometry()[i_node].Coordinates()[i_dim] -= delta;
                }else{
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = 0.0;
                }
            }
        }

        //Added by inigo

        const int wake = pGetPrimalElement()->GetValue(WAKE);

        if (wake == 0) // Normal element (non-wake) - eventually an embedded
        {
            BoundedMatrix<double, NumNodes, Dim> x;
            for(unsigned int i_node = 0; i_node < NumNodes; i_node++){
                    x( i_node , 0 ) = GetGeometry()[i_node].X();
                    x( i_node , 1 ) = GetGeometry()[i_node].Y();
            }
            auto p = PotentialFlow::GetPotentialOnNormalElement<2,3>(*pGetPrimalElement());
            BoundedMatrix<double, Dim*NumNodes, NumNodes> test = ZeroMatrix(Dim*NumNodes, NumNodes);

            const double crOutput0 =             x(0,0) - x(1,0);
            const double crOutput1 =             -x(2,1);
            const double crOutput2 =             crOutput1 + x(0,1);
            const double crOutput3 =             -x(2,0);
            const double crOutput4 =             crOutput3 + x(0,0);
            const double crOutput5 =             x(0,1) - x(1,1);
            const double crOutput6 =             crOutput0*crOutput2 - crOutput4*crOutput5;
            const double crOutput7 =             pow(crOutput6, -2);
            const double crOutput8 =             crOutput3 + x(1,0);
            const double crOutput9 =             -p[2];
            const double crOutput10 =             crOutput6*(crOutput9 + p[1]);
            const double crOutput11 =             crOutput1 + x(1,1);
            const double crOutput12 =             crOutput0*crOutput8 + crOutput11*crOutput5;
            const double crOutput13 =             crOutput11*crOutput2 + crOutput4*crOutput8;
            const double crOutput14 =             crOutput12*p[2] - crOutput13*p[1] + p[0]*(pow(crOutput11, 2) + pow(crOutput8, 2));
            const double crOutput15 =             crOutput8*p[0];
            const double crOutput16 =             crOutput4*p[1];
            const double crOutput17 =             2*crOutput16;
            const double crOutput18 =             -2*x(0,0) + x(1,0) + x(2,0);
            const double crOutput19 =             crOutput0*crOutput4 + crOutput2*crOutput5;
            const double crOutput20 =             crOutput13*p[0] + crOutput19*p[2] - p[1]*(pow(crOutput2, 2) + pow(crOutput4, 2));
            const double crOutput21 =             crOutput0*p[2];
            const double crOutput22 =             2*crOutput21;
            const double crOutput23 =             crOutput12*p[0] - crOutput19*p[1] + p[2]*(pow(crOutput0, 2) + pow(crOutput5, 2));
            const double crOutput24 =             crOutput11*p[0];
            const double crOutput25 =             crOutput2*p[1];
            const double crOutput26 =             2*crOutput25;
            const double crOutput27 =             -2*x(0,1) + x(1,1) + x(2,1);
            const double crOutput28 =             crOutput5*p[2];
            const double crOutput29 =             2*crOutput28;
            const double crOutput30 =             2*crOutput15;
            const double crOutput31 =             x(0,0) - 2*x(1,0) + x(2,0);
            const double crOutput32 =             crOutput6*(crOutput9 + p[0]);
            const double crOutput33 =             2*crOutput24;
            const double crOutput34 =             x(0,1) - 2*x(1,1) + x(2,1);
            const double crOutput35 =             x(0,0) + x(1,0) - 2*x(2,0);
            const double crOutput36 =             crOutput6*(p[0] - p[1]);
            const double crOutput37 =             x(0,1) + x(1,1) - 2*x(2,1);
            test(0,0)=crOutput7*(crOutput10*crOutput8 + crOutput11*crOutput14);
            test(0,1)=-crOutput7*(crOutput11*crOutput20 + crOutput6*(-crOutput15 + crOutput17 + crOutput18*p[2]));
            test(0,2)=crOutput7*(crOutput11*crOutput23 - crOutput6*(crOutput15 + crOutput18*p[1] + crOutput22));
            test(1,0)=crOutput7*(crOutput10*crOutput11 - crOutput14*crOutput8);
            test(1,1)=crOutput7*(crOutput20*crOutput8 - crOutput6*(-crOutput24 + crOutput26 + crOutput27*p[2]));
            test(1,2)=-crOutput7*(crOutput23*crOutput8 + crOutput6*(crOutput24 + crOutput27*p[1] + crOutput29));
            test(2,0)=-crOutput7*(crOutput14*crOutput2 + crOutput6*(-crOutput16 + crOutput30 + crOutput31*p[2]));
            test(2,1)=crOutput7*(crOutput2*crOutput20 + crOutput32*crOutput4);
            test(2,2)=-crOutput7*(crOutput2*crOutput23 + crOutput6*(crOutput16 - crOutput22 + crOutput31*p[0]));
            test(3,0)=crOutput7*(crOutput14*crOutput4 - crOutput6*(-crOutput25 + crOutput33 + crOutput34*p[2]));
            test(3,1)=crOutput7*(crOutput2*crOutput32 - crOutput20*crOutput4);
            test(3,2)=crOutput7*(crOutput23*crOutput4 - crOutput6*(crOutput25 - crOutput29 + crOutput34*p[0]));
            test(4,0)=crOutput7*(crOutput14*crOutput5 + crOutput6*(crOutput21 + crOutput30 - crOutput35*p[1]));
            test(4,1)=-crOutput7*(crOutput20*crOutput5 + crOutput6*(-crOutput17 + crOutput21 + crOutput35*p[0]));
            test(4,2)=crOutput7*(crOutput0*crOutput36 + crOutput23*crOutput5);
            test(5,0)=-crOutput7*(crOutput0*crOutput14 - crOutput6*(crOutput28 + crOutput33 - crOutput37*p[1]));
            test(5,1)=crOutput7*(crOutput0*crOutput20 - crOutput6*(-crOutput26 + crOutput28 + crOutput37*p[0]));
            test(5,2)=crOutput7*(-crOutput0*crOutput23 + crOutput36*crOutput5);

            //array_1d<double, NumNodes> phis = pGetPrimalElement()->GetPotentialOnNormalElement(pGetPrimalElement());
            

            if(pGetPrimalElement()->Id()==100)
            {
            KRATOS_WATCH(rOutput(1,1))
            KRATOS_WATCH(test(1,1))
            KRATOS_WATCH(p)
            KRATOS_WATCH(delta)
            }
        }


        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) 
    {
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0)//normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                    else
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                }
            }
        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }
        }


    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) 
    {
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }

    template <class TPrimalElement>
    int AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) 
    {

        KRATOS_TRY

        int Check = mpPrimalElement -> Check(rCurrentProcessInfo);

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
            {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,
                                                    this->GetGeometry()[i]);

                return Check;
            }
        }

        return 0;

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    template <class TPrimalElement>
    std::string AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Info() const 
    {
        std::stringstream buffer;
        buffer << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << Id();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::PrintData(std::ostream& rOStream) const 
    {
        pGetGeometry()->PrintData(rOStream);
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    /*PROTECTED*/

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    /*PRIVATE*/

    template <class TPrimalElement>
    double AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetPerturbationSize()
    {
        const double delta = this->GetValue(SCALE_FACTOR);
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const 
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer) 
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointAnalyticalIncompressiblePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    //template class AdjointAnalyticalIncompressiblePotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.

