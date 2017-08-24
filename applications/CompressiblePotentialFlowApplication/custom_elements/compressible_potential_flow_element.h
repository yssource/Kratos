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

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED

// #define SYMMETRIC_CONSTRAINT_APPLICATION

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"
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
class CompressiblePotentialFlowElement : public Element
{
public:

    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double,TNumNodes> phis, distances;
        double rho;
        double vol;

        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
    };



    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(CompressiblePotentialFlowElement);

    typedef Condition::WeakPointer ConditionWeakPointerType;
    
    typedef Condition::Pointer ConditionPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    CompressiblePotentialFlowElement(IndexType NewId = 0) {};

    /**
     * Constructor using an array of nodes
     */
    CompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId, ThisNodes) {};

    /**
     * Constructor using Geometry
     */
    CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId, pGeometry) {};

    /**
     * Constructor using Properties
     */
    CompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Element(NewId, pGeometry, pProperties) {};

    /**
     * Copy Constructor
     */
    CompressiblePotentialFlowElement(CompressiblePotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~CompressiblePotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CompressiblePotentialFlowElement & operator=(CompressiblePotentialFlowElement const& rOther)
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, pGeom, pProperties));
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
        return Element::Pointer(new CompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
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
                rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();

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
                    rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
                else
                    rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
                else
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE,0).EquationId();
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
                rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
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
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
                else
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
                else
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
            }
        }
    }

    //get Area*normal from the condition
    void Calculate(const Variable< array_1d<double,3> >& rVariable,
        array_1d<double,3>& Output,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        ConditionPointerType pCond = pGetCondition();//the search of conditions needs to be implemented
        pCond->Calculate(VELOCITY,Output, rCurrentProcessInfo);
    }

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector& values, int Step = 0) override
    {
        //gather nodal data
        for(unsigned int i=0; i<NumNodes; i++)
        {
            values[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        }
    }


    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        ElementalData<NumNodes,Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        //gather nodal data
        for(unsigned int i=0; i<NumNodes; i++)
        {
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        }
        
        
        //TEST:
        bool kutta_element = false;
        for(unsigned int i=0; i<NumNodes; ++i)
            if(GetGeometry()[i].Is(STRUCTURE))
            {
                kutta_element = true;
                break;
            }

        if(this->IsNot(MARKER))//normal element (non-wake) - eventually an embedded
        {
            if(rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
                rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);
            if(rRightHandSideVector.size() != NumNodes)
                rRightHandSideVector.resize(NumNodes,false);
            rLeftHandSideMatrix.clear();

            ComputeLHSGaussPointContribution(data.vol,rLeftHandSideMatrix,data);
                        
            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
            
        }
        else //it is a wake element
        {
            GetWakeDistances(data.distances);
            
            //note that the lhs and rhs have double the size!!
            if(rLeftHandSideMatrix.size1() != 2*NumNodes || rLeftHandSideMatrix.size2() != 2*NumNodes)
                rLeftHandSideMatrix.resize(2*NumNodes,2*NumNodes,false);
            if(rRightHandSideVector.size() != 2*NumNodes)
                rRightHandSideVector.resize(2*NumNodes,false);
            rLeftHandSideMatrix.clear();
            
            //subdivide the element
            constexpr unsigned int nvolumes = 3*(Dim-1);
            bounded_matrix<double,NumNodes, Dim > Points;
            array_1d<double,nvolumes> Volumes;
            bounded_matrix<double, nvolumes, NumNodes > GPShapeFunctionValues;
            array_1d<double,nvolumes> PartitionsSign;
            std::vector<Matrix> GradientsValue(nvolumes);
            bounded_matrix<double,nvolumes, 2> NEnriched;
            
            for(unsigned int i=0; i<GradientsValue.size(); ++i)
                GradientsValue[i].resize(2,Dim,false);
           
            
            
            for(unsigned int i = 0; i<NumNodes; ++i)
            {
                const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
                for(unsigned int k = 0; k<Dim; ++k)
                {
                    Points(i, k) = coords[k];
                }
            }
            
            const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(Points,
                                                                                            data.DN_DX,
                                                                                            data.distances,
                                                                                            Volumes, 
                                                                                            GPShapeFunctionValues, 
                                                                                            PartitionsSign, 
                                                                                            GradientsValue, 
                                                                                            NEnriched);
            //compute the lhs and rhs that would correspond to it not being divided
            Matrix lhs_positive = ZeroMatrix(NumNodes,NumNodes);
            Matrix lhs_negative = ZeroMatrix(NumNodes,NumNodes);
            
            for(unsigned int i=0; i<nsubdivisions; ++i)
            {
                if(PartitionsSign[i] > 0)
                    ComputeLHSGaussPointContribution(Volumes[i],lhs_positive,data);
                else
                    ComputeLHSGaussPointContribution(Volumes[i],lhs_negative,data);
            }

            //also next version works - NON SYMMETRIC - but it does not require a penalty
//                 array_1d<double,Dim> n = prod(data.DN_DX,data.distances); //rCurrentProcessInfo[VELOCITY]; 
//                 n /= norm_2(n);
//                 bounded_matrix<double,Dim,Dim> nn = outer_prod(n,n);
//                 bounded_matrix<double,NumNodes,Dim> tmp = prod(data.DN_DX,nn);
//                 bounded_matrix<double,NumNodes,NumNodes> constraint = data.vol*prod(tmp, trans(data.DN_DX));
//                                 
//                 bounded_matrix<double,Dim,Dim> P = IdentityMatrix(Dim,Dim) - nn;
//                 noalias(tmp) = prod(data.DN_DX,P);
//                 bounded_matrix<double,NumNodes,NumNodes> tangent_constraint = /*1e3**/data.vol*prod(tmp, trans(data.DN_DX));
                if(kutta_element == true)
                {
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        for(unsigned int j=0; j<NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i,j)            =  lhs_positive(i,j); 
                            rLeftHandSideMatrix(i,j+NumNodes)   =  0.0; 
                            
                            rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  lhs_negative(i,j); 
                            rLeftHandSideMatrix(i+NumNodes,j)          =  0.0; 
                        }
                    }
                }
                else
                {
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        for(unsigned int j=0; j<NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i,j)                   =  lhs_positive(i,j); 
                            rLeftHandSideMatrix(i,j+NumNodes)          =  0.0; 
                            
                            rLeftHandSideMatrix(i+NumNodes,j+NumNodes) =  lhs_negative(i,j); 
                            rLeftHandSideMatrix(i+NumNodes,j)          =  0.0; 
                        }
                    }
                    
                
                    //side1  -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {
                        if(data.distances[i]<0)
                        {
                                for(unsigned int j=0; j<NumNodes; ++j)
                                {
                                    rLeftHandSideMatrix(i,j)          = lhs_positive(i,j); 
                                    rLeftHandSideMatrix(i,j+NumNodes) = -lhs_positive(i,j); 
                                }
                        }
                    }
                    
                    //side2 -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs
                    for(unsigned int i=0; i<NumNodes; ++i)
                    {                            
                        if(data.distances[i]>0)
                        {   
                            for(unsigned int j=0; j<NumNodes; ++j)
                                {
                                    rLeftHandSideMatrix(i+NumNodes,j+NumNodes) = lhs_negative(i,j);
                                    rLeftHandSideMatrix(i+NumNodes,j) = -lhs_negative(i,j); 
                                }
                        }
                    }
                }
            Vector split_element_values(NumNodes*2);
            GetValuesOnSplitElement(split_element_values, data.distances);
            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,split_element_values);
        }
        
    }


    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        //TODO: improve speed
        Matrix tmp;
        CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
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
            KRATOS_THROW_ERROR(std::logic_error, "CompressiblePotentialFlowElement found with Id 0 or negative","")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on CompressiblePotentialFlowElement -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0","")
        }

        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE ) == false )
                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable POSITIVE_FACE_PRESSURE on node ", this->GetGeometry()[i].Id() )
            }


        return 0;

        KRATOS_CATCH("");
    }

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if(rValues.size() != 1) rValues.resize(1);

        if (rVariable == PRESSURE)
        {
            double p = 0.0;

            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            if(active && !this->Is(MARKER))
            {
                const array_1d<double,3> vinfinity = rCurrentProcessInfo[VELOCITY];
                const double vinfinity_norm2 = inner_prod(vinfinity,vinfinity);

                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                //gather nodal data
                for(unsigned int i=0; i<NumNodes; i++)
                {
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                }

                const array_1d<double,Dim> v = prod(trans(data.DN_DX), data.phis);


                p = (vinfinity_norm2 - inner_prod(v,v))/vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));
            }


            rValues[0] = p;
        }
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if(rValues.size() != 1) rValues.resize(1);

        if (rVariable == VELOCITY)
        {
            bool active = true;
            if ((this)->IsDefined(ACTIVE))
                active = (this)->Is(ACTIVE);

            array_1d<double,3> v = ZeroVector();
            if(this->IsNot(MARKER) && active==true)
            {
                ElementalData<NumNodes,Dim> data;

                //calculate shape functions
                GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

                //gather nodal data
                for(unsigned int i=0; i<NumNodes; i++)
                {
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                }

                array_1d<double,Dim> vaux = -prod(trans(data.DN_DX), data.phis);
                
                for(unsigned int k=0; k<Dim; k++) v[k] = vaux[k];
            }



            rValues[0] = v;
        }
    }

    /**
     * @brief Calculates the adjoint matrix for potential.
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * potential transposed
     */    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        //TODO: improve speed
        Vector tmp;
        CalculateLocalSystem(rLeftHandSideMatrix, tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(-rLeftHandSideMatrix); // transpose
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
            const unsigned int dim = GetGeometry().WorkingSpaceDimension();
            const unsigned int nnodes = GetGeometry().size();

            //matrix of coordinates
            bounded_matrix<double,NumNodes, Dim> x(NumNodes,dim);
            for(unsigned int i=0; i<NumNodes; ++i)
                for(unsigned int k=0; k<dim; k++)
                    x(i,k) = GetGeometry()[i].Coordinates()[k];
            
            
            bounded_matrix<double, NumNodes, Dim > DN;   
            bounded_matrix<double,NumNodes, NumNodes*Dim> DRDx;
            //std::cout << "DIM #" << Dim;
            array_1d<double,NumNodes> N;
            double vol;
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN, N, vol);
            
            //gather nodal data
            array_1d<double,NumNodes> p;
            for(unsigned int i=0; i<nnodes; i++)
                p[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);

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
                
                rOutput = trans(DRDx);

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
        buffer << "CompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressiblePotentialFlowElement #" << Id();
    }

/// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
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
    void GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
//         for(unsigned int i = 0; i<NumNodes; i++)
//             distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }
    
    
    void ComputeLHSGaussPointContribution(
        const double weight,
        Matrix& lhs,
        const ElementalData<NumNodes,Dim>& data)
    {
        noalias(lhs) += weight*prod(data.DN_DX, trans(data.DN_DX));
    }

    void ComputeRHSGaussPointContribution(
        const double weight,
        Vector& rhs,
        const ElementalData<NumNodes,Dim>& data)
    {
        array_1d<double,Dim> grad = prod(trans(data.DN_DX), data.phis);
        noalias(rhs) -= weight*prod(data.DN_DX, grad);
    }


    void GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
        }
    }




    ///@}
    ///@name Protected Operations
    ///@{

    ConditionPointerType pGetCondition()
	{
	   return mpCondition.lock();
	}

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
    ConditionWeakPointerType mpCondition;


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
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

}; // Class CompressiblePotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
