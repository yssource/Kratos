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
    virtual ~CompressiblePotentialFlowElement() {};

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
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
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
    virtual void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
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

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    virtual void GetValuesVector(Vector& values, int Step = 0)
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
    virtual void CalculateLocalSystem(
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
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override
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

    virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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

    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
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

    
    virtual void Calculate(const Variable<Matrix >& rVariable,
            Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int nnodes = GetGeometry().size();
        if(rVariable == ADJOINT_MATRIX_1)
        {
            Matrix lhs;
            Vector rhs;
            CalculateLocalSystem(lhs,rhs,const_cast<ProcessInfo&>(rCurrentProcessInfo)); //TODO: the const cast here is HORRIBLE
            rOutput = trans(lhs);
        }
        else if(rVariable == SHAPE_DERIVATIVE_MATRIX_1)
        {
            //matrix of coordinates
            bounded_matrix<double,Dim, NumNodes> x(dim,NumNodes);
            for(unsigned int i=0; i<NumNodes; ++i)
                for(unsigned int k=0; k<dim; k++)
                    x(k,i) = GetGeometry()[i].Coordinates()[k];
            
            bounded_matrix<double, NumNodes, Dim > DN;   
            bounded_matrix<double,NumNodes, Dim> DRDx;
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
                
                
                const double cDRDx0 =             DN(0,1)*x(0,0) + DN(1,1)*x(1,0) + DN(2,1)*x(2,0);
                const double cDRDx1 =             DN(0,0)*cDRDx0;
                const double cDRDx2 =             DN(0,0)*x(0,0) + DN(1,0)*x(1,0) + DN(2,0)*x(2,0);
                const double cDRDx3 =             -DN(0,1)*cDRDx2 + cDRDx1;
                const double cDRDx4 =             DN(0,1)*x(0,1) + DN(1,1)*x(1,1) + DN(2,1)*x(2,1);
                const double cDRDx5 =             DN(0,0)*x(0,1) + DN(1,0)*x(1,1) + DN(2,0)*x(2,1);
                const double cDRDx6 =             DN(0,1)*cDRDx5;
                const double cDRDx7 =             DN(0,0)*cDRDx4 - cDRDx6;
                const double cDRDx8 =             cDRDx0*cDRDx5;
                const double cDRDx9 =             cDRDx2*cDRDx4 - cDRDx8;
                const double cDRDx10 =             cDRDx3*p[0];
                const double cDRDx11 =             DN(1,0)*cDRDx0;
                const double cDRDx12 =             -DN(1,1)*cDRDx2 + cDRDx11;
                const double cDRDx13 =             cDRDx12*p[1];
                const double cDRDx14 =             DN(2,0)*cDRDx0;
                const double cDRDx15 =             -DN(2,1)*cDRDx2 + cDRDx14;
                const double cDRDx16 =             cDRDx15*p[2];
                const double cDRDx17 =             cDRDx10 + cDRDx13 + cDRDx16;
                const double cDRDx18 =             cDRDx17/pow(cDRDx9, 3);
                const double cDRDx19 =             pow(cDRDx9, -2);
                const double cDRDx20 =             cDRDx19*cDRDx3;
                const double cDRDx21 =             1.0/cDRDx9;
                const double cDRDx22 =             DN(0,0)*DN(1,1);
                const double cDRDx23 =             DN(0,1)*DN(1,0);
                const double cDRDx24 =             DN(1,0)*cDRDx0*cDRDx21;
                const double cDRDx25 =             DN(1,1)*cDRDx2*cDRDx21;
                const double cDRDx26 =             cDRDx22 - cDRDx23 + cDRDx24*cDRDx7 - cDRDx25*cDRDx7;
                const double cDRDx27 =             DN(0,0)*DN(2,1);
                const double cDRDx28 =             DN(0,1)*DN(2,0);
                const double cDRDx29 =             DN(2,0)*cDRDx0*cDRDx21;
                const double cDRDx30 =             DN(2,1)*cDRDx2*cDRDx21;
                const double cDRDx31 =             cDRDx27 - cDRDx28 + cDRDx29*cDRDx7 - cDRDx30*cDRDx7;
                const double cDRDx32 =             cDRDx10*cDRDx21*cDRDx7 + cDRDx26*p[1] + cDRDx31*p[2];
                const double cDRDx33 =             DN(0,1)*cDRDx19*cDRDx5;
                const double cDRDx34 =             1.0/cDRDx2;
                const double cDRDx35 =             DN(0,0)*cDRDx34;
                const double cDRDx36 =             cDRDx21*cDRDx6;
                const double cDRDx37 =             cDRDx21*cDRDx34*cDRDx5;
                const double cDRDx38 =             cDRDx0*cDRDx19*cDRDx5;
                const double cDRDx39 =             -cDRDx1*cDRDx37 - cDRDx35 + cDRDx36 - cDRDx38*cDRDx7;
                const double cDRDx40 =             cDRDx33*cDRDx7 + cDRDx35*cDRDx39;
                const double cDRDx41 =             cDRDx21*cDRDx8 + 1;
                const double cDRDx42 =             cDRDx35*cDRDx41 - cDRDx36;
                const double cDRDx43 =             DN(1,1)*cDRDx5;
                const double cDRDx44 =             cDRDx21*cDRDx43;
                const double cDRDx45 =             DN(1,0)*cDRDx34;
                const double cDRDx46 =             cDRDx41*cDRDx45 - cDRDx44;
                const double cDRDx47 =             DN(2,1)*cDRDx5;
                const double cDRDx48 =             cDRDx21*cDRDx47;
                const double cDRDx49 =             DN(2,0)*cDRDx34;
                const double cDRDx50 =             cDRDx41*cDRDx49 - cDRDx48;
                const double cDRDx51 =             cDRDx42*p[0] + cDRDx46*p[1] + cDRDx50*p[2];
                const double cDRDx52 =             DN(1,1)*cDRDx19*cDRDx5;
                const double cDRDx53 =             cDRDx39*cDRDx45 + cDRDx52*cDRDx7;
                const double cDRDx54 =             DN(2,1)*cDRDx19*cDRDx5;
                const double cDRDx55 =             cDRDx39*cDRDx49 + cDRDx54*cDRDx7;
                const double cDRDx56 =             cDRDx40*p[0] + cDRDx53*p[1] + cDRDx55*p[2];
                const double cDRDx57 =             cDRDx17*cDRDx19;
                const double cDRDx58 =             cDRDx12*cDRDx19;
                const double cDRDx59 =             cDRDx15*cDRDx19;
                const double cDRDx60 =             cDRDx21*rho;
                const double cDRDx61 =             2*cDRDx17*cDRDx19;
                const double cDRDx62 =             DN(0,0)*cDRDx0*cDRDx34;
                const double cDRDx63 =             cDRDx21*cDRDx5;
                const double cDRDx64 =             DN(0,0) + cDRDx3*cDRDx63;
                const double cDRDx65 =             DN(0,0)*DN(0,1) + cDRDx3*cDRDx36 - cDRDx62*cDRDx64;
                const double cDRDx66 =             DN(1,0)*cDRDx0*cDRDx34;
                const double cDRDx67 =             cDRDx22 + cDRDx3*cDRDx44 - cDRDx64*cDRDx66;
                const double cDRDx68 =             DN(2,0)*cDRDx0*cDRDx34;
                const double cDRDx69 =             cDRDx27 + cDRDx3*cDRDx48 - cDRDx64*cDRDx68;
                const double cDRDx70 =             cDRDx65*p[0] + cDRDx67*p[1] + cDRDx69*p[2];
                const double cDRDx71 =             2*cDRDx17*cDRDx19*cDRDx3;
                const double cDRDx72 =             -cDRDx12*cDRDx71;
                const double cDRDx73 =             -cDRDx15*cDRDx71;
                const double cDRDx74 =             DN(1,0)*cDRDx4 - cDRDx43;
                const double cDRDx75 =             DN(0,0)*cDRDx0*cDRDx21;
                const double cDRDx76 =             DN(0,1)*cDRDx2*cDRDx21;
                const double cDRDx77 =             -cDRDx22 + cDRDx23 + cDRDx74*cDRDx75 - cDRDx74*cDRDx76;
                const double cDRDx78 =             DN(1,0)*DN(2,1);
                const double cDRDx79 =             DN(1,1)*DN(2,0);
                const double cDRDx80 =             cDRDx29*cDRDx74 - cDRDx30*cDRDx74 + cDRDx78 - cDRDx79;
                const double cDRDx81 =             cDRDx13*cDRDx21*cDRDx74 + cDRDx77*p[0] + cDRDx80*p[2];
                const double cDRDx82 =             -cDRDx11*cDRDx37 - cDRDx38*cDRDx74 + cDRDx44 - cDRDx45;
                const double cDRDx83 =             cDRDx33*cDRDx74 + cDRDx35*cDRDx82;
                const double cDRDx84 =             cDRDx45*cDRDx82 + cDRDx52*cDRDx74;
                const double cDRDx85 =             cDRDx49*cDRDx82 + cDRDx54*cDRDx74;
                const double cDRDx86 =             cDRDx83*p[0] + cDRDx84*p[1] + cDRDx85*p[2];
                const double cDRDx87 =             DN(1,0) + cDRDx12*cDRDx63;
                const double cDRDx88 =             cDRDx12*cDRDx36 + cDRDx23 - cDRDx62*cDRDx87;
                const double cDRDx89 =             DN(1,0)*DN(1,1) + cDRDx12*cDRDx44 - cDRDx66*cDRDx87;
                const double cDRDx90 =             cDRDx12*cDRDx48 - cDRDx68*cDRDx87 + cDRDx78;
                const double cDRDx91 =             cDRDx88*p[0] + cDRDx89*p[1] + cDRDx90*p[2];
                const double cDRDx92 =             -cDRDx12*cDRDx15*cDRDx61;
                const double cDRDx93 =             DN(2,0)*cDRDx4 - cDRDx47;
                const double cDRDx94 =             -cDRDx27 + cDRDx28 + cDRDx75*cDRDx93 - cDRDx76*cDRDx93;
                const double cDRDx95 =             cDRDx24*cDRDx93 - cDRDx25*cDRDx93 - cDRDx78 + cDRDx79;
                const double cDRDx96 =             cDRDx16*cDRDx21*cDRDx93 + cDRDx94*p[0] + cDRDx95*p[1];
                const double cDRDx97 =             -cDRDx14*cDRDx37 - cDRDx38*cDRDx93 + cDRDx48 - cDRDx49;
                const double cDRDx98 =             cDRDx33*cDRDx93 + cDRDx35*cDRDx97;
                const double cDRDx99 =             cDRDx45*cDRDx97 + cDRDx52*cDRDx93;
                const double cDRDx100 =             cDRDx49*cDRDx97 + cDRDx54*cDRDx93;
                const double cDRDx101 =             cDRDx100*p[2] + cDRDx98*p[0] + cDRDx99*p[1];
                const double cDRDx102 =             DN(2,0) + cDRDx15*cDRDx63;
                const double cDRDx103 =             -cDRDx102*cDRDx62 + cDRDx15*cDRDx36 + cDRDx28;
                const double cDRDx104 =             -cDRDx102*cDRDx66 + cDRDx15*cDRDx44 + cDRDx79;
                const double cDRDx105 =             DN(2,0)*DN(2,1) - cDRDx102*cDRDx68 + cDRDx15*cDRDx48;
                const double cDRDx106 =             cDRDx103*p[0] + cDRDx104*p[1] + cDRDx105*p[2];
                DRDx(0,0)=rho*(cDRDx18*cDRDx3*cDRDx7 + cDRDx20*cDRDx32 - cDRDx40*cDRDx51 - cDRDx42*cDRDx56);
                DRDx(0,1)=rho*(cDRDx26*cDRDx57 + cDRDx32*cDRDx58 - cDRDx46*cDRDx56 - cDRDx51*cDRDx53);
                DRDx(0,2)=rho*(cDRDx31*cDRDx57 + cDRDx32*cDRDx59 - cDRDx50*cDRDx56 - cDRDx51*cDRDx55);
                DRDx(1,0)=cDRDx60*(-pow(cDRDx3, 2)*cDRDx61 + cDRDx42*cDRDx70 + cDRDx51*cDRDx65);
                DRDx(1,1)=cDRDx60*(cDRDx46*cDRDx70 + cDRDx51*cDRDx67 + cDRDx72);
                DRDx(1,2)=cDRDx60*(cDRDx50*cDRDx70 + cDRDx51*cDRDx69 + cDRDx73);
                DRDx(2,0)=rho*(cDRDx20*cDRDx81 - cDRDx42*cDRDx86 - cDRDx51*cDRDx83 + cDRDx57*cDRDx77);
                DRDx(2,1)=rho*(cDRDx12*cDRDx18*cDRDx74 - cDRDx46*cDRDx86 - cDRDx51*cDRDx84 + cDRDx58*cDRDx81);
                DRDx(2,2)=rho*(-cDRDx50*cDRDx86 - cDRDx51*cDRDx85 + cDRDx57*cDRDx80 + cDRDx59*cDRDx81);
                DRDx(3,0)=cDRDx60*(cDRDx42*cDRDx91 + cDRDx51*cDRDx88 + cDRDx72);
                DRDx(3,1)=cDRDx60*(-pow(cDRDx12, 2)*cDRDx61 + cDRDx46*cDRDx91 + cDRDx51*cDRDx89);
                DRDx(3,2)=cDRDx60*(cDRDx50*cDRDx91 + cDRDx51*cDRDx90 + cDRDx92);
                DRDx(4,0)=rho*(-cDRDx101*cDRDx42 + cDRDx20*cDRDx96 - cDRDx51*cDRDx98 + cDRDx57*cDRDx94);
                DRDx(4,1)=rho*(-cDRDx101*cDRDx46 - cDRDx51*cDRDx99 + cDRDx57*cDRDx95 + cDRDx58*cDRDx96);
                DRDx(4,2)=rho*(-cDRDx100*cDRDx51 - cDRDx101*cDRDx50 + cDRDx15*cDRDx18*cDRDx93 + cDRDx59*cDRDx96);
                DRDx(5,0)=cDRDx60*(cDRDx103*cDRDx51 + cDRDx106*cDRDx42 + cDRDx73);
                DRDx(5,1)=cDRDx60*(cDRDx104*cDRDx51 + cDRDx106*cDRDx46 + cDRDx92);
                DRDx(5,2)=cDRDx60*(cDRDx105*cDRDx51 + cDRDx106*cDRDx50 - pow(cDRDx15, 2)*cDRDx61);

                
                
                rOutput = trans(DRDx);

            }
            else
            {
                KRATOS_ERROR << "sorry, SHAPE_DERIVATIVE_MATRIX not yet implemented in 3D";
            }
        }
        
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "CompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CompressiblePotentialFlowElement #" << Id();
    }

/// Print object's data.

    void PrintData(std::ostream& rOStream) const
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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer)
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
