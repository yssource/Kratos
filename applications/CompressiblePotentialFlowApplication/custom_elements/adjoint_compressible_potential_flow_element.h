//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_ADJOINT_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_ADJOINT_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED

// #define SYMMETRIC_CONSTRAINT_APPLICATION

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
//#include "applications/CompressiblePotentialFlowApplication/custom_elements/compressible_potential_flow_element.h" 
#include "custom_elements/compressible_potential_flow_element.h" 

namespace Kratos
{

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

template< int Dim, int NumNodes >
class AdjointCompressiblePotentialFlowElement : public CompressiblePotentialFlowElement< Dim, NumNodes>
{
public:

    
    ///@name Type Definitions
    ///@{

    typedef CompressiblePotentialFlowElement< Dim, NumNodes> BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointCompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(AdjointCompressiblePotentialFlowElement);

    typedef typename BaseType::ConditionWeakPointerType ConditionWeakPointerType;
    
    typedef typename BaseType::ConditionPointerType ConditionPointerType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::EquationIdVectorType EquationIdVectorType;

    typedef Element::DofsVectorType DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    AdjointCompressiblePotentialFlowElement(IndexType NewId = 0) {};

    /**
     * Constructor using an array of nodes
     */
    AdjointCompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes):CompressiblePotentialFlowElement< Dim, NumNodes>(NewId, ThisNodes) {};

    /**
     * Constructor using Geometry
     */
    AdjointCompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry):CompressiblePotentialFlowElement< Dim, NumNodes>(NewId, pGeometry) {};

    /**
     * Constructor using Properties
     */
    AdjointCompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):CompressiblePotentialFlowElement< Dim, NumNodes>(NewId, pGeometry, pProperties) {};

    /**
     * Copy Constructor
     */
    AdjointCompressiblePotentialFlowElement(AdjointCompressiblePotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~AdjointCompressiblePotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdjointCompressiblePotentialFlowElement & operator=(AdjointCompressiblePotentialFlowElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointCompressiblePotentialFlowElement(NewId, Element::GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointCompressiblePotentialFlowElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointCompressiblePotentialFlowElement(NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties()));
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            for (unsigned int i = 0; i < NumNodes; i++)
                rResult[i] = Element::GetGeometry()[i].GetDof(ADJOINT_POSITIVE_FACE_PRESSURE).EquationId();

        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            BaseType::GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = Element::GetGeometry()[i].GetDof(ADJOINT_POSITIVE_FACE_PRESSURE).EquationId();
                else
                    rResult[i] = Element::GetGeometry()[i].GetDof(ADJOINT_NEGATIVE_FACE_PRESSURE).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = Element::GetGeometry()[i].GetDof(ADJOINT_POSITIVE_FACE_PRESSURE).EquationId();
                else
                    rResult[NumNodes+i] = Element::GetGeometry()[i].GetDof(ADJOINT_NEGATIVE_FACE_PRESSURE).EquationId();
            }
        }


    }


    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            for (unsigned int i = 0; i < NumNodes; i++)
                rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(ADJOINT_POSITIVE_FACE_PRESSURE);
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances;
            BaseType::GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(ADJOINT_POSITIVE_FACE_PRESSURE);
                else
                    rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(ADJOINT_NEGATIVE_FACE_PRESSURE);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = Element::GetGeometry()[i].pGetDof(ADJOINT_POSITIVE_FACE_PRESSURE);
                else
                    rElementalDofList[NumNodes+i] = Element::GetGeometry()[i].pGetDof(ADJOINT_NEGATIVE_FACE_PRESSURE);
            }
        }
    }

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
     void GetValuesVector(Vector& values, int Step = 0) override
     {
         if(this->IsNot(MARKER))//normal element
         {
            if(values.size()!=NumNodes)
               values.resize(NumNodes,false);
            
            //gather nodal data
            for(unsigned int i=0; i<NumNodes; i++)
            {
                values[i] = Element::GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_POSITIVE_FACE_PRESSURE);//careful with the sizes!
            }
         }
         else//wake element
         {
            if(values.size()!=2*NumNodes)
            {
                values.resize(2*NumNodes,false);
            }

            array_1d<double,NumNodes> distances;
            BaseType::GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    values[i] = Element::GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_POSITIVE_FACE_PRESSURE);
                else
                    values[i] = Element::GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_NEGATIVE_FACE_PRESSURE);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    values[NumNodes+i] = Element::GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_POSITIVE_FACE_PRESSURE);
                else
                    values[NumNodes+i] = Element::GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_NEGATIVE_FACE_PRESSURE);
            }

         }
         
     }

   
    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "AdjointCompressiblePotentialFlowElement found with Id 0 or negative","")
        }

        if (Element::GetGeometry().Area() <= 0)
        {
            std::cout << "error on AdjointCompressiblePotentialFlowElement -> " << Element::Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0","")
        }

        for ( unsigned int i = 0; i < Element::GetGeometry().size(); i++ )
        {
            if ( Element::GetGeometry()[i].SolutionStepsDataHas( ADJOINT_POSITIVE_FACE_PRESSURE ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable ADJOINT_POSITIVE_FACE_PRESSURE on node ", Element::GetGeometry()[i].Id() )
            }


        return BaseType::Check(rCurrentProcessInfo);

        KRATOS_CATCH("");
    }



    /**
     * Calculate the transposed gradient of the element's residual w.r.t. design variable (i.e. in this case w.r.t. potential).
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            const unsigned int dim = Element::GetGeometry().WorkingSpaceDimension();
            const unsigned int nnodes = Element::GetGeometry().size();

            //matrix of coordinates
            bounded_matrix<double,NumNodes, Dim> x(NumNodes,dim);
            for(unsigned int i=0; i<NumNodes; ++i)
                for(unsigned int k=0; k<dim; k++)
                    x(i,k) = Element::GetGeometry()[i].Coordinates()[k];
            
            
            bounded_matrix<double, NumNodes, Dim > DN;   
            
            //std::cout << "DIM #" << Dim;
            array_1d<double,NumNodes> N;
            double vol;
            GeometryUtils::CalculateGeometryData(Element::GetGeometry(), DN, N, vol);
            
            //gather nodal data
            array_1d<double,NumNodes> p;
            for(unsigned int i=0; i<nnodes; i++)
                p[i] = Element::GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);

            const double rho  = 1; //clearly this should be changed...
            
            if(dim == 2)
            {
                const double cDRDx0 =             x(0,0) - x(1,0);
                const double cDRDx1 =             -x(2,1);
                const double cDRDx2 =             cDRDx1 + x(0,1);
                const double cDRDx3 =             -x(2,0);
                const double cDRDx4 =             cDRDx3 + x(0,0);
                const double cDRDx5 =             x(0,1) - x(1,1);
                const double cDRDx6 =             cDRDx0*cDRDx2 - cDRDx4*cDRDx5;
                const double cDRDx7 =             0.5*rho/pow(cDRDx6, 2);
                const double cDRDx8 =             cDRDx3 + x(1,0);
                const double cDRDx9 =             -p[2];
                const double cDRDx10 =             cDRDx9 + p[1];
                const double cDRDx11 =             cDRDx10*cDRDx6;
                const double cDRDx12 =             cDRDx1 + x(1,1);
                const double cDRDx13 =             cDRDx0*p[2] - cDRDx4*p[1] + cDRDx8*p[0];
                const double cDRDx14 =             cDRDx12*p[0] - cDRDx2*p[1] + cDRDx5*p[2];
                const double cDRDx15 =             cDRDx12*cDRDx14 + cDRDx13*cDRDx8;
                const double cDRDx16 =             cDRDx9 + p[0];
                const double cDRDx17 =             p[0] - p[1];
                const double cDRDx18 =             cDRDx13*cDRDx4 + cDRDx14*cDRDx2;
                const double cDRDx19 =             cDRDx16*cDRDx6;
                const double cDRDx20 =             cDRDx0*cDRDx13 + cDRDx14*cDRDx5;
                const double cDRDx21 =             cDRDx17*cDRDx6;
                if(this->IsNot(MARKER))//normal element
                {
                    bounded_matrix<double,NumNodes, NumNodes*Dim> DRDx;

                    DRDx(0,0)=cDRDx7*(cDRDx11*cDRDx8 + cDRDx12*cDRDx15);
                    DRDx(0,1)=cDRDx7*(cDRDx11*cDRDx12 - cDRDx15*cDRDx8);
                    DRDx(0,2)=-cDRDx7*(cDRDx15*cDRDx2 + cDRDx6*(cDRDx13 + cDRDx16*cDRDx8));
                    DRDx(0,3)=cDRDx7*(cDRDx15*cDRDx4 - cDRDx6*(cDRDx12*cDRDx16 + cDRDx14));
                    DRDx(0,4)=cDRDx7*(cDRDx15*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx8));
                    DRDx(0,5)=-cDRDx7*(cDRDx0*cDRDx15 - cDRDx6*(cDRDx12*cDRDx17 + cDRDx14));
                    DRDx(1,0)=-cDRDx7*(cDRDx12*cDRDx18 - cDRDx6*(-cDRDx10*cDRDx4 + cDRDx13));
                    DRDx(1,1)=cDRDx7*(cDRDx18*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx2 + cDRDx14));
                    DRDx(1,2)=cDRDx7*(cDRDx18*cDRDx2 + cDRDx19*cDRDx4);
                    DRDx(1,3)=cDRDx7*(-cDRDx18*cDRDx4 + cDRDx19*cDRDx2);
                    DRDx(1,4)=-cDRDx7*(cDRDx18*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx4));
                    DRDx(1,5)=cDRDx7*(cDRDx0*cDRDx18 - cDRDx6*(cDRDx14 + cDRDx17*cDRDx2));
                    DRDx(2,0)=cDRDx7*(cDRDx12*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx10 + cDRDx13));
                    DRDx(2,1)=-cDRDx7*(cDRDx20*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx5 + cDRDx14));
                    DRDx(2,2)=-cDRDx7*(cDRDx2*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx16 + cDRDx13));
                    DRDx(2,3)=cDRDx7*(cDRDx20*cDRDx4 + cDRDx6*(cDRDx14 - cDRDx16*cDRDx5));
                    DRDx(2,4)=cDRDx7*(cDRDx0*cDRDx21 + cDRDx20*cDRDx5);
                    DRDx(2,5)=cDRDx7*(-cDRDx0*cDRDx20 + cDRDx21*cDRDx5);

                    rOutput = trans(DRDx);
                }
                else//wake element
                {
                    bounded_matrix<double,2*NumNodes, NumNodes*Dim> DRDx;

                    DRDx(0,0)=cDRDx7*(cDRDx11*cDRDx8 + cDRDx12*cDRDx15);
                    DRDx(0,1)=cDRDx7*(cDRDx11*cDRDx12 - cDRDx15*cDRDx8);
                    DRDx(0,2)=-cDRDx7*(cDRDx15*cDRDx2 + cDRDx6*(cDRDx13 + cDRDx16*cDRDx8));
                    DRDx(0,3)=cDRDx7*(cDRDx15*cDRDx4 - cDRDx6*(cDRDx12*cDRDx16 + cDRDx14));
                    DRDx(0,4)=cDRDx7*(cDRDx15*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx8));
                    DRDx(0,5)=-cDRDx7*(cDRDx0*cDRDx15 - cDRDx6*(cDRDx12*cDRDx17 + cDRDx14));
                    DRDx(1,0)=-cDRDx7*(cDRDx12*cDRDx18 - cDRDx6*(-cDRDx10*cDRDx4 + cDRDx13));
                    DRDx(1,1)=cDRDx7*(cDRDx18*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx2 + cDRDx14));
                    DRDx(1,2)=cDRDx7*(cDRDx18*cDRDx2 + cDRDx19*cDRDx4);
                    DRDx(1,3)=cDRDx7*(-cDRDx18*cDRDx4 + cDRDx19*cDRDx2);
                    DRDx(1,4)=-cDRDx7*(cDRDx18*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx4));
                    DRDx(1,5)=cDRDx7*(cDRDx0*cDRDx18 - cDRDx6*(cDRDx14 + cDRDx17*cDRDx2));
                    DRDx(2,0)=cDRDx7*(cDRDx12*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx10 + cDRDx13));
                    DRDx(2,1)=-cDRDx7*(cDRDx20*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx5 + cDRDx14));
                    DRDx(2,2)=-cDRDx7*(cDRDx2*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx16 + cDRDx13));
                    DRDx(2,3)=cDRDx7*(cDRDx20*cDRDx4 + cDRDx6*(cDRDx14 - cDRDx16*cDRDx5));
                    DRDx(2,4)=cDRDx7*(cDRDx0*cDRDx21 + cDRDx20*cDRDx5);
                    DRDx(2,5)=cDRDx7*(-cDRDx0*cDRDx20 + cDRDx21*cDRDx5);

                    //For wake elements simply a repetition is implemented CHECK!!
                    DRDx(3,0)=cDRDx7*(cDRDx11*cDRDx8 + cDRDx12*cDRDx15);
                    DRDx(3,1)=cDRDx7*(cDRDx11*cDRDx12 - cDRDx15*cDRDx8);
                    DRDx(3,2)=-cDRDx7*(cDRDx15*cDRDx2 + cDRDx6*(cDRDx13 + cDRDx16*cDRDx8));
                    DRDx(3,3)=cDRDx7*(cDRDx15*cDRDx4 - cDRDx6*(cDRDx12*cDRDx16 + cDRDx14));
                    DRDx(3,4)=cDRDx7*(cDRDx15*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx8));
                    DRDx(3,5)=-cDRDx7*(cDRDx0*cDRDx15 - cDRDx6*(cDRDx12*cDRDx17 + cDRDx14));
                    DRDx(4,0)=-cDRDx7*(cDRDx12*cDRDx18 - cDRDx6*(-cDRDx10*cDRDx4 + cDRDx13));
                    DRDx(4,1)=cDRDx7*(cDRDx18*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx2 + cDRDx14));
                    DRDx(4,2)=cDRDx7*(cDRDx18*cDRDx2 + cDRDx19*cDRDx4);
                    DRDx(4,3)=cDRDx7*(-cDRDx18*cDRDx4 + cDRDx19*cDRDx2);
                    DRDx(4,4)=-cDRDx7*(cDRDx18*cDRDx5 + cDRDx6*(cDRDx13 + cDRDx17*cDRDx4));
                    DRDx(4,5)=cDRDx7*(cDRDx0*cDRDx18 - cDRDx6*(cDRDx14 + cDRDx17*cDRDx2));
                    DRDx(5,0)=cDRDx7*(cDRDx12*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx10 + cDRDx13));
                    DRDx(5,1)=-cDRDx7*(cDRDx20*cDRDx8 + cDRDx6*(-cDRDx10*cDRDx5 + cDRDx14));
                    DRDx(5,2)=-cDRDx7*(cDRDx2*cDRDx20 - cDRDx6*(-cDRDx0*cDRDx16 + cDRDx13));
                    DRDx(5,3)=cDRDx7*(cDRDx20*cDRDx4 + cDRDx6*(cDRDx14 - cDRDx16*cDRDx5));
                    DRDx(5,4)=cDRDx7*(cDRDx0*cDRDx21 + cDRDx20*cDRDx5);
                    DRDx(5,5)=cDRDx7*(-cDRDx0*cDRDx20 + cDRDx21*cDRDx5);

                    rOutput = trans(DRDx);
                }


                
                
                // const double cDRDx0 =             DN(0,0)*x(0,0) + DN(1,0)*x(1,0) + DN(2,0)*x(2,0);
                // const double cDRDx1 =             DN(0,1)*x(0,1) + DN(1,1)*x(1,1) + DN(2,1)*x(2,1);
                // const double cDRDx2 =             DN(0,0)*x(0,1) + DN(1,0)*x(1,1) + DN(2,0)*x(2,1);
                // const double cDRDx3 =             DN(0,1)*x(0,0) + DN(1,1)*x(1,0) + DN(2,1)*x(2,0);
                // const double cDRDx4 =             cDRDx0*cDRDx1 - cDRDx2*cDRDx3;
                // const double cDRDx5 =             0.5*rho/pow(cDRDx4, 2);
                // const double cDRDx6 =             DN(0,0)*cDRDx3 - DN(0,1)*cDRDx0;
                // const double cDRDx7 =             DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0);
                // const double cDRDx8 =             DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0);
                // const double cDRDx9 =             cDRDx7*p[1] + cDRDx8*p[2];
                // const double cDRDx10 =             cDRDx4*cDRDx9;
                // const double cDRDx11 =             DN(0,0)*cDRDx1 - DN(0,1)*cDRDx2;
                // const double cDRDx12 =             DN(1,0)*cDRDx3 - DN(1,1)*cDRDx0;
                // const double cDRDx13 =             DN(2,0)*cDRDx3 - DN(2,1)*cDRDx0;
                // const double cDRDx14 =             cDRDx12*p[1] + cDRDx13*p[2] + cDRDx6*p[0];
                // const double cDRDx15 =             DN(1,0)*cDRDx1 - DN(1,1)*cDRDx2;
                // const double cDRDx16 =             DN(2,0)*cDRDx1 - DN(2,1)*cDRDx2;
                // const double cDRDx17 =             cDRDx11*p[0] + cDRDx15*p[1] + cDRDx16*p[2];
                // const double cDRDx18 =             cDRDx11*cDRDx17 + cDRDx14*cDRDx6;
                // const double cDRDx19 =             DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0);
                // const double cDRDx20 =             -cDRDx19*p[2] + cDRDx7*p[0];
                // const double cDRDx21 =             cDRDx14*cDRDx7;
                // const double cDRDx22 =             cDRDx17*cDRDx7;
                // const double cDRDx23 =             cDRDx19*p[1] + cDRDx8*p[0];
                // const double cDRDx24 =             cDRDx14*cDRDx8;
                // const double cDRDx25 =             cDRDx17*cDRDx8;
                // const double cDRDx26 =             cDRDx12*cDRDx14 + cDRDx15*cDRDx17;
                // const double cDRDx27 =             cDRDx20*cDRDx4;
                // const double cDRDx28 =             cDRDx14*cDRDx19;
                // const double cDRDx29 =             cDRDx17*cDRDx19;
                // const double cDRDx30 =             cDRDx13*cDRDx14 + cDRDx16*cDRDx17;
                // const double cDRDx31 =             cDRDx23*cDRDx4;
                // DRDx(0,0)=cDRDx5*(cDRDx10*cDRDx6 + cDRDx11*cDRDx18);
                // DRDx(0,1)=-cDRDx5*(-cDRDx10*cDRDx11 + cDRDx18*cDRDx6);
                // DRDx(0,2)=cDRDx5*(cDRDx15*cDRDx18 - cDRDx4*(cDRDx20*cDRDx6 + cDRDx21));
                // DRDx(0,3)=-cDRDx5*(cDRDx12*cDRDx18 + cDRDx4*(cDRDx11*cDRDx20 + cDRDx22));
                // DRDx(0,4)=cDRDx5*(cDRDx16*cDRDx18 - cDRDx4*(cDRDx23*cDRDx6 + cDRDx24));
                // DRDx(0,5)=-cDRDx5*(cDRDx13*cDRDx18 + cDRDx4*(cDRDx11*cDRDx23 + cDRDx25));
                // DRDx(1,0)=cDRDx5*(cDRDx11*cDRDx26 + cDRDx4*(cDRDx12*cDRDx9 + cDRDx21));
                // DRDx(1,1)=-cDRDx5*(cDRDx26*cDRDx6 - cDRDx4*(cDRDx15*cDRDx9 + cDRDx22));
                // DRDx(1,2)=-cDRDx5*(cDRDx12*cDRDx27 - cDRDx15*cDRDx26);
                // DRDx(1,3)=-cDRDx5*(cDRDx12*cDRDx26 + cDRDx15*cDRDx27);
                // DRDx(1,4)=cDRDx5*(cDRDx16*cDRDx26 - cDRDx4*(cDRDx12*cDRDx23 + cDRDx28));
                // DRDx(1,5)=-cDRDx5*(cDRDx13*cDRDx26 + cDRDx4*(cDRDx15*cDRDx23 + cDRDx29));
                // DRDx(2,0)=cDRDx5*(cDRDx11*cDRDx30 + cDRDx4*(cDRDx13*cDRDx9 + cDRDx24));
                // DRDx(2,1)=-cDRDx5*(cDRDx30*cDRDx6 - cDRDx4*(cDRDx16*cDRDx9 + cDRDx25));
                // DRDx(2,2)=cDRDx5*(cDRDx15*cDRDx30 + cDRDx4*(-cDRDx13*cDRDx20 + cDRDx28));
                // DRDx(2,3)=-cDRDx5*(cDRDx12*cDRDx30 - cDRDx4*(-cDRDx16*cDRDx20 + cDRDx29));
                // DRDx(2,4)=-cDRDx5*(cDRDx13*cDRDx31 - cDRDx16*cDRDx30);
                // DRDx(2,5)=-cDRDx5*(cDRDx13*cDRDx30 + cDRDx16*cDRDx31);
                
                // const double cDRDx0 =             DN(0,0)*x(0,0) + DN(1,0)*x(1,0) + DN(2,0)*x(2,0);
                // const double cDRDx1 =             DN(0,1)*x(0,1) + DN(1,1)*x(1,1) + DN(2,1)*x(2,1);
                // const double cDRDx2 =             DN(0,0)*x(0,1) + DN(1,0)*x(1,1) + DN(2,0)*x(2,1);
                // const double cDRDx3 =             DN(0,1)*x(0,0) + DN(1,1)*x(1,0) + DN(2,1)*x(2,0);
                // const double cDRDx4 =             1.0/(cDRDx0*cDRDx1 - cDRDx2*cDRDx3);
                // const double cDRDx5 =             0.5*cDRDx4*rho;
                // const double cDRDx6 =             DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0);
                // const double cDRDx7 =             DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0);
                // const double cDRDx8 =             cDRDx6*p[1] + cDRDx7*p[2];
                // const double cDRDx9 =             DN(0,0)*cDRDx3 - DN(0,1)*cDRDx0;
                // const double cDRDx10 =             DN(0,0)*cDRDx1 - DN(0,1)*cDRDx2;
                // const double cDRDx11 =             cDRDx10*cDRDx4;
                // const double cDRDx12 =             DN(1,0)*cDRDx3 - DN(1,1)*cDRDx0;
                // const double cDRDx13 =             DN(2,0)*cDRDx3 - DN(2,1)*cDRDx0;
                // const double cDRDx14 =             cDRDx12*p[1] + cDRDx13*p[2] + cDRDx9*p[0];
                // const double cDRDx15 =             DN(1,0)*cDRDx1 - DN(1,1)*cDRDx2;
                // const double cDRDx16 =             DN(2,0)*cDRDx1 - DN(2,1)*cDRDx2;
                // const double cDRDx17 =             cDRDx10*p[0] + cDRDx15*p[1] + cDRDx16*p[2];
                // const double cDRDx18 =             cDRDx10*cDRDx17 + cDRDx14*cDRDx9;
                // const double cDRDx19 =             cDRDx4*cDRDx9;
                // const double cDRDx20 =             DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0);
                // const double cDRDx21 =             -cDRDx20*p[2] + cDRDx6*p[0];
                // const double cDRDx22 =             cDRDx14*cDRDx6;
                // const double cDRDx23 =             cDRDx15*cDRDx4;
                // const double cDRDx24 =             cDRDx17*cDRDx6;
                // const double cDRDx25 =             cDRDx12*cDRDx4;
                // const double cDRDx26 =             cDRDx20*p[1] + cDRDx7*p[0];
                // const double cDRDx27 =             cDRDx14*cDRDx7;
                // const double cDRDx28 =             cDRDx16*cDRDx4;
                // const double cDRDx29 =             cDRDx17*cDRDx7;
                // const double cDRDx30 =             cDRDx13*cDRDx4;
                // const double cDRDx31 =             cDRDx12*cDRDx14 + cDRDx15*cDRDx17;
                // const double cDRDx32 =             cDRDx14*cDRDx20;
                // const double cDRDx33 =             cDRDx17*cDRDx20;
                // const double cDRDx34 =             cDRDx13*cDRDx14 + cDRDx16*cDRDx17;
                // DRDx(0,0)=cDRDx5*(cDRDx11*cDRDx18 + cDRDx8*cDRDx9);
                // DRDx(0,1)=cDRDx5*(cDRDx10*cDRDx8 - cDRDx18*cDRDx19);
                // DRDx(0,2)=cDRDx5*(cDRDx18*cDRDx23 - cDRDx21*cDRDx9 - cDRDx22);
                // DRDx(0,3)=-cDRDx5*(cDRDx10*cDRDx21 + cDRDx18*cDRDx25 + cDRDx24);
                // DRDx(0,4)=cDRDx5*(cDRDx18*cDRDx28 - cDRDx26*cDRDx9 - cDRDx27);
                // DRDx(0,5)=-cDRDx5*(cDRDx10*cDRDx26 + cDRDx18*cDRDx30 + cDRDx29);
                // DRDx(1,0)=cDRDx5*(cDRDx11*cDRDx31 + cDRDx12*cDRDx8 + cDRDx22);
                // DRDx(1,1)=cDRDx5*(cDRDx15*cDRDx8 - cDRDx19*cDRDx31 + cDRDx24);
                // DRDx(1,2)=cDRDx5*(-cDRDx12*cDRDx21 + cDRDx23*cDRDx31);
                // DRDx(1,3)=-cDRDx5*(cDRDx15*cDRDx21 + cDRDx25*cDRDx31);
                // DRDx(1,4)=cDRDx5*(-cDRDx12*cDRDx26 + cDRDx28*cDRDx31 - cDRDx32);
                // DRDx(1,5)=-cDRDx5*(cDRDx15*cDRDx26 + cDRDx30*cDRDx31 + cDRDx33);
                // DRDx(2,0)=cDRDx5*(cDRDx11*cDRDx34 + cDRDx13*cDRDx8 + cDRDx27);
                // DRDx(2,1)=cDRDx5*(cDRDx16*cDRDx8 - cDRDx19*cDRDx34 + cDRDx29);
                // DRDx(2,2)=cDRDx5*(-cDRDx13*cDRDx21 + cDRDx23*cDRDx34 + cDRDx32);
                // DRDx(2,3)=-cDRDx5*(cDRDx16*cDRDx21 + cDRDx25*cDRDx34 - cDRDx33);
                // DRDx(2,4)=cDRDx5*(-cDRDx13*cDRDx26 + cDRDx28*cDRDx34);
                // DRDx(2,5)=-cDRDx5*(cDRDx16*cDRDx26 + cDRDx30*cDRDx34);              
                
                

            }
            else
            {
                KRATOS_ERROR << "sorry, SHAPE_SENSITIVITY not yet implemented in 3D";
            }            
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
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
        std::stringstream buffer;
        buffer << "AdjointCompressiblePotentialFlowElement #" << Element::Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AdjointCompressiblePotentialFlowElement #" << Element::Id();
    }

/// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        Element::pGetGeometry()->PrintData(rOStream);
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


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointCompressiblePotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_ADJOINT_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
