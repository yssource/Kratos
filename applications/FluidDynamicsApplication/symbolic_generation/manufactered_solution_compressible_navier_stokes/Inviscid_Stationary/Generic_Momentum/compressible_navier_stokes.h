//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES) 
#define  KRATOS_COMPRESSIBLE_NAVIER_STOKES

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h" 

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

// TODO: UPDATE THIS INFORMATION
/**this element is a 3D stokes element, stabilized by employing an ASGS stabilization
* formulation is described in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwZ2Zxd09YUmlPZ28/view?usp=sharing
* symbolic implementation is defined in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwaXRKRUpDbmx4VXM/view?usp=sharing
*/
template< unsigned int TDim, unsigned int BlockSize = TDim+2, unsigned int TNumNodes = TDim + 1 >
class CompressibleNavierStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokes);
    struct ElementDataStruct
    {
        bounded_matrix<double, TNumNodes, BlockSize> U, Un, Unn;
//         array_1d<double, TNumNodes > c;            //speed of sound
//         array_1d<double, TNumNodes > tau1,tau2,tau3; // stabilization matrix parameters
        bounded_matrix<double, TNumNodes, TDim> f_ext;
        array_1d<double,TNumNodes> r; // At the moment considering all parameters as constant in the domain (mu, nu, etc...)
        array_1d<double, TDim> f_gauss;
        double r_gauss;

        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        
        double bdf0;
        double bdf1;
        double bdf2;
        double h;             // Element size 
        double volume;        // In 2D: element area. In 3D: element volume
        double mu;
        double nu;
        double lambda;          
        double cv;
        double y;               //gamma
        double c;               // TO DO : temporarily use for testing
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    CompressibleNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    CompressibleNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~CompressibleNavierStokes() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< CompressibleNavierStokes < TDim,BlockSize, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< CompressibleNavierStokes < TDim,BlockSize, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        constexpr unsigned int MatrixSize = TNumNodes*(BlockSize);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);
        
        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Gauss point position
        //bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        //GetShapeFunctionsOnGauss(Ncontainer);
        unsigned int ngauss = GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_4).size1();
        
        for(unsigned int igauss = 0; igauss<ngauss; igauss++)
        {
            //noalias(data.N) = row(Ncontainer, igauss);
            noalias(data.N) = row(GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_4), igauss);
            double xgauss = 0;
            double ygauss = 0;
            
            for(unsigned int i = 0; i<GetGeometry().size(); i++){
                xgauss += GetGeometry()[i].X()*data.N[i];
                ygauss += GetGeometry()[i].Y()*data.N[i];
            }
            double x = xgauss;
            double y = ygauss;
            /*
            //mu = 0.3 lambda = 1.0
            data.f_gauss[0] =  -0.111111111111111*(132.75*pow(x,7) + 1858.5*pow(x,6)*y + 5403.6*pow(x,6) + 10690.2*pow(x,5)*y*y + 56986.56*pow(x,5)*y + 49872.24*pow(x,5) + 32562.0*pow(x,4)*pow(y,3)\
             + 241963.2*pow(x,4)*y*y + 435854.88*pow(x,4)*y + 210706.2*pow(x,4) + 55908.0*pow(x,3)*pow(y,4) + 520819.2*pow(x,3)*pow(y,3) + 1480152.96*pow(x,3)*y*y + 1484427.456*pow(x,3)*y + 484003.584*pow(x,3)\
              + 52344.0*x*x*pow(y,5) + 579859.2*x*x*pow(y,4) + 2410179.84*x*x*pow(y,3) + 3830741.568*x*x*y*y + 2581992.3456*x*x*y + 627573.10464*x*x + 22608.0*x*pow(y,6) + 290995.2*x*pow(y,5) + 1836460.8*x*pow(y,4) + 4252435.2*x*pow(y,3) \
              + 4504826.88*x*y*y + 2252610.10944*x*y + 433458.118656*x + 2246.39999999999*pow(y,7) + 35435.52*pow(y,6) + 495659.52*pow(y,5) + 1686024.576*pow(y,4) + 2553713.0496*pow(y,3) + 1990095.962112*y*y + 784439.0240256*y + 124394.8253184)/(pow((0.5*x + 1.0*y + 0.8),6)*pow((1.0*x + 2.0*y + 1.6),2));

            data.f_gauss[1] =  -0.111111111111111*(689.85*pow(x,7) + 5971.5*pow(x,6)*y + 7634.16*pow(x,6) + 17397.0*pow(x,5)*y*y + 52488.0*pow(x,5)*y + 33500.16*pow(x,5) + 8837.99999999999*pow(x,4)*pow(y,3)\
             + 96321.6*pow(x,4)*y*y + 163447.2*pow(x,4)*y + 73767.492*pow(x,4) - 56052.0000000001*pow(x,3)*pow(y,4)  - 107481.6*pow(x,3)*pow(y,3) + 61943.0400000002*pow(x,3)*y*y + 191723.04*pow(x,3)*y + 81663.2064*pow(x,3) - 126244.8*x*x*pow(y,5)\
              - 589766.4*x*x*pow(y,4) - 871061.759999999*x*x*pow(y,3) - 469086.624*x*x*y*y - 27392.2560000003*x*x*y + 32518.434816*x*x - 103824.0*x*pow(y,6) - 720276.48*x*pow(y,5) - 1677265.92*x*pow(y,4) - 1816463.232*x*pow(y,3)\
               - 968731.5456*x*y*y - 226401.73056*x*y - 12957.6075264*x - 29664.0*pow(y,7) - 291456.0*pow(y,6) - 917890.56*pow(y,5)- 1403075.52*pow(y,4) - 1174588.416*pow(y,3) - 544222.49472*y*y - 127797.16608*y - 11216.8009728)/(pow((0.5*x + 1.0*y + 0.8),6)*pow((1.0*x + 2.0*y + 1.6),2));


            data.r_gauss = 1.0*(1876765122.0*pow(x,10)+ 23089082702.4*pow(x,9)*y - 63938302295.4*pow(x,9) + 128351247537.6*pow(x,8)*y*y - 1172895707099.52*pow(x,8)*y - 1228800908084.49*pow(x,8) + 484876304409.6*pow(x,7)*pow(y,3) - 8876629673389.44*pow(x,7)*y*y - 18740620748962.4*pow(x,7)*y - 8804661442182.97*pow(x,7) \
            + 1622861303616.0*pow(x,6)*pow(y,4) - 36385870333478.4*pow(x,6)*pow(y,3) - 120723116005349.0*pow(x,6)*y*y - 116351930845578.0*pow(x,6)*y - 35114886793212.3*pow(x,6)+ 4832728185600.0*pow(x,5)*pow(y,5)\
             - 87551210790086.4*pow(x,5)*pow(y,4) - 427420063774781.0*pow(x,5)*pow(y,3) - 643138151885957.0*pow(x,5)*y*y - 397598988705811.0*pow(x,5)*y - 87865282098847.8*pow(x,5) + 10839029848012.8*pow(x,4)*pow(y,6)\
              - 121972239344794.0*pow(x,4)*pow(y,5) - 901930470081722.0*pow(x,4)*pow(y,4) - 1.92210955940701e+15*pow(x,4)*pow(y,3) - 1.84085360786283e+15*pow(x,4)*y*y - 830959475264704.0*pow(x,4)*y - 144119895079452.0*pow(x,4)\
               + 16235828722483.2*pow(x,3)*pow(y,7) - 82652913395619.8*pow(x,3)*pow(y,6) - 1.14160337838682e+15*pow(x,3)*pow(y,5) - 3.33704892748786e+15*pow(x,3)*pow(y,4) - 4.45171310153756e+15*pow(x,3)*pow(y,3) - 3.09531007610314e+15*pow(x,3)*y*y - 1.09377456693309e+15*pow(x,3)*y - 155451801761574.0*pow(x,3) + 15015397825843.2*x*x*pow(y,8)\
                + 2034550356664.3*x*x*pow(y,7) - 815032165073662.0*x*x*pow(y,67) - 3.33498537624359e+15*x*x*pow(y,5) - 5.91013828661146e+15*x*x*pow(y,4) - 5.66809904212905e+15*x*x*pow(y,3) - 3.07228737453247e+15*x*x*y*y - 887782996821446.0*x*x*y - 106481250928374.0*x*x + 7670601197568.0*x*pow(y,9) + 37870000416921.6*x*pow(y,8)\
                - 270390367655784.0*x*pow(y,7) - 1.74689052152098e+15*x*pow(y,6) - 4.06145054743746e+15*x*pow(y,5) - 5.09015589584126e+15*x*pow(y,4) - 3.78104783025142e+15*x*pow(y,3) - 1.67073533535682e+15*x*y*y - 406758187091081.0*x*y - 42052786842382.1*x + 1631992264704.0*pow(y,10) + 16307569867161.6*pow(y,9) - 17980620652149.1*pow(y,8) - 357371958291277.0*pow(y,7)\
                 - 1.11820991874469e+15*pow(y,6) - 1.78673365491928e+15*pow(y,5) - 1.71722437656645e+15*pow(y,4) - 1.03517743872642e+15*pow(y,3) - 384480096619006.0*y*y - 80581883473861.1*y - 7295005034606.26)/(9061.9453125*pow(x,11) + 199362.796875*pow(x,10)*y + 159490.2375*pow(x,10) + 1993627.96875*pow(x,9)*y*y + 3189804.75*pow(x,9)*y +\
                  1275921.9*pow(x,9) + 11961767.8125*pow(x,8)*pow(y,3) + 28708242.75*pow(x,8)*y*y + 22966594.2*pow(x,8)*y + 6124425.12*pow(x,8) + 47847071.25*pow(x,7)*pow(y,4) + 153110628.0*pow(x,7)*pow(y,3) + 183732753.6*pow(x,7)*y*y + 97990801.92*pow(x,7)*y + 19598160.384*pow(x,7) + 133971799.5*pow(x,6)*pow(y,5) + 535887198.0*pow(x,6)*pow(x,7)*pow(y,4) + 857419516.8*pow(x,6)*pow(y,3) + 685935613.44*pow(x,6)*y*y + 274374245.376*pow(x,6)*y +\
                   43899879.26016*pow(x,6)+ 267943599.0*pow(x,5)*pow(y,6) + 1286129275.2*pow(x,5)*pow(y,5) + 2572258550.4*pow(x,5)*pow(y,4) + 2743742453.76*pow(x,5)*pow(y,3) + 1646245472.256*pow(x,5)*y*y + 526798551.12192*pow(x,5)*y + 70239806.816256*pow(x,5) + 382776570.0*pow(x,4)*pow(y,7) + 2143548792.0*pow(x,4)*pow(y,6) + 5144517100.8*pow(x,4)*pow(y,5) + 6859356134.4*pow(x,4)*pow(y,4) + 5487484907.52*pow(x,4)*pow(y,3) + 2633992755.6096*pow(x,4)*y*y + 702398068.16256*pow(x,4)*y + 80274064.932864*pow(x,4) + \
                   382776570.0*pow(x,3)*pow(y,8) + 2449770048.0*pow(x,3)*pow(y,7) + 6859356134.4*pow(x,3)*pow(y,6) + 10974969815.04*pow(x,3)*pow(y,5) + 10974969815.04*pow(x,3)*pow(y,4) + 7023980681.6256*pow(x,3)*pow(y,3) + 2809592272.65024*pow(x,3)*y*y + 642192519.462912*pow(x,3)*y + 64219251.9462912*pow(x,3) + 255184380.0*x*x*pow(y,9) + 1837327536.0*x*x*pow(y,8) + 5879448115.2*x*x*pow(y,7) + 10974969815.04*x*x*pow(y,6) + 13169963778.048*x*x*pow(y,5) + 10535971022.4384*x*x*pow(y,4)\
                    + 5619184545.30048*x*x*pow(y,3) + 1926577558.38874*x*x*y*y+ 385315511.677747*x*x*y + 34250267.7046887*x*x + 102073752.0*x*pow(y,10) + 816590016.0*x*pow(y,9) + \
                   2939724057.6*x*pow(y,8) + 6271411322.88*x*pow(y,7) + 8779975852.032*x*pow(y,6) + 8428776817.95072*x*pow(y,5) + 5619184545.30048*x*pow(y,4) + 2568770077.85165*x*pow(y,3) + 770631023.355495*x*y*y + 137001070.818755*x*y \
                   + 10960085.6655004*x + 18558864.0*pow(y,11) + 163318003.2*pow(y,10) + 653272012.8*pow(y,9) + 1567852830.72*pow(y,8) + 2508564529.152*pow(y,7) + 2809592272.65024*pow(y,6) + 2247673818.12019*pow(y,5) + 1284385038.92582*pow(y,4) + 513754015.57033*pow(y,3) + 137001070.818755*y*y + 21920171.3310007*y + 1594194.27861824);
            */
            
            
            //mu = 0.0 lambda = 0.0
            data.f_gauss[0] = -1.0*(118.0*pow(x,4)+ 944.0*pow(x,3)*y + 4236.8*pow(x,3) + 2422.4*x*x*y*y + 18437.12*x*x*y + 23084.8*x*x + 2137.6*x*pow(y,3) + 21596.16*x*y*y + 69991.424*x*y + 43266.048*x + 249.599999999999*pow(y,4) + 3338.24000000001*pow(y,3) + 46595.072*y*y + 68653.8752*y + 26709.1968)/(pow((0.5*x + 1.0*y + 0.8),3)*pow((1.0*x + 2.0*y + 1.6),2));

            data.f_gauss[1] = -1.0*(613.2*pow(x,4) + 1628.8*pow(x,3)*y + 3842.56*pow(x,3) - 1667.2*x*x*y*y + 4008.96000000001*x*x*y + 6573.056*x*x - 6592.0*x*pow(y,3) - 19589.12*x*y*y - 9490.43199999999*x*y + 1569.5872*x - 3296.0*pow(y,4) - 24473.6*pow(y,3) - 36884.48*y*y - 18776.064*y - 2595.2256)/(pow((0.5*x + 1.0*y + 0.8),3)*pow((1.0*x + 2.0*y + 1.6),2));
            
            data.r_gauss = 1.0*(12944.0*pow(x,5) + 29804.8*pow(x,4)*y - 544532.8*pow(x,4) + 69427.2*pow(x,3)*y*y - 3710957.44*pow(x,3)*y - 4449297.28*pow(x,3) + 422195.2*x*x*pow(y,3) - 7279169.28*x*x*y*y - 22964259.072*x*x*y - 11712234.496*x*x + 773888.0*x*pow(y,4) - 3016870.4*x*pow(y,3) - 35006937.6*x*y*y - 41496424.448*x*y - 12654321.664*x + 351744.0*pow(y,5) + 2107801.6*pow(y,4) - 14556523.52*pow(y,3)\
                 - 34091323.392*y*y - 23001186.304*y - 4783538.176)/(0.0625*pow(x,6) + 0.75*pow(x,5)*y + 0.6*pow(x,5) + 3.75*pow(x,4)*y*y + 6.0*pow(x,4)*y + 2.4*pow(x,4) + 10.0*pow(x,3)*pow(y,3) + 24.0*pow(x,3)*y*y + 19.2*pow(x,3)*y + 5.12*pow(x,3) + 15.0*x*x*pow(y,4) + 48.0*x*x*pow(y,3) + 57.6*x*x*y*y + 30.72*x*x*y + 6.144*x*x + 12.0*x*pow(y,5) + 48.0*x*pow(y,4) + 76.8*x*pow(y,3) + 61.44*x*y*y + 24.576*x*y + 3.93216*x \
                 + 4.0*pow(y,6) + 19.2*pow(y,5) + 38.4*pow(y,4) + 40.96*pow(y,3) + 24.576*y*y + 7.86432*y + 1.048576);
            
            

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            // here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix  *= data.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        std::cout<<this->Id()<<" "<<rRightHandSideVector<<std::endl;
        for(unsigned int i=0; i<GetGeometry().size(); i++){
            std::cout<<GetGeometry()[i].Id()<<" "<<GetGeometry()[i].Coordinates()<<std::endl;
            std::cout<<GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM)<<std::endl;
            std::cout<<GetGeometry()[i].FastGetSolutionStepValue(DENSITY)<<std::endl;
            std::cout<<GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY)<<std::endl;
        }

        KRATOS_CATCH("Error in Compressible Navier Stokes Element Symbolic")
    }


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(BlockSize);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);
 
        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Gauss point position
        //bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        //GetShapeFunctionsOnGauss(Ncontainer);

        unsigned int ngauss = GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_4).size1();

        // Loop on gauss point
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<ngauss; igauss++)
        {
            //noalias(data.N) = row(Ncontainer, igauss);
            noalias(data.N) = row(GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_4), igauss);


            ComputeGaussPointRHSContribution(rhs_local, data);
            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rRightHandSideVector) += rhs_local;
        }

        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("")

    }

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        if(MOMENTUM.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MOMENTUM Key is 0. Check if the application was correctly registered.","");
        if(TOTAL_ENERGY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"TOTAL_ENERGY Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(DYNAMIC_VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.","");
        if(KINEMATIC_VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"KINEMATIC_VISCOSITY Key is 0. Check if the application was correctly registered.","");
        if(CONDUCTIVITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONDUCTIVITY Key is 0. Check if the application was correctly registered.","");
        if(SPECIFIC_HEAT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"SPECIFIC_HEAT Key is 0. Check if the application was correctly registered.","");
        if(HEAT_CAPACITY_RATIO.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"HEAT_CAPACITY_RATIO Key is 0. Check if the application was correctly registered.","");
        
    

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing MOMENTUM variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TOTAL_ENERGY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) == false ||
               this->GetGeometry()[i].HasDofFor(MOMENTUM_Y) == false ||
               this->GetGeometry()[i].HasDofFor(MOMENTUM_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing MOMENTUM component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TOTAL_ENERGY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(DENSITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing DENSITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

       return 0;
	
        KRATOS_CATCH("");
    }
    

    void Calculate(const Variable<double>& rVariable,
                           double& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        if (rVariable == ERROR_RATIO)
        {
//             rOutput = this->SubscaleErrorEstimate(data);
            this->SetValue(ERROR_RATIO, rOutput);
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
        return "CompressibleNavierStokes #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const override

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static member variables
    ///@{
    

    ///@}
    ///@name Protected member Variables
    ///@{

    // Symbolic function implementing the element
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*(BlockSize),TNumNodes*(BlockSize)>& lhs, const ElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(BlockSize)>& rhs, const ElementDataStruct& data);

    double SubscaleErrorEstimate(const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    CompressibleNavierStokes() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    // Auxiliar function to fill the element data structure
    void FillElementData(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {   
        // Getting data for the given geometry
        // double Volume; // In 2D cases Volume variable contains the element area
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), rData.DN_DX, rData.N, rData.volume);

        // Compute element size
        rData.h = ComputeH(rData.DN_DX);

        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];
        
        Properties& r_properties = this->GetProperties();
        rData.nu = r_properties.GetValue(KINEMATIC_VISCOSITY);
        rData.mu =  r_properties.GetValue(DYNAMIC_VISCOSITY);
        rData.lambda = r_properties.GetValue(CONDUCTIVITY);
        rData.cv = r_properties.GetValue(SPECIFIC_HEAT);
        rData.y = r_properties.GetValue(HEAT_CAPACITY_RATIO);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            //const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& moment = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& moment_n = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM,1);
            const array_1d<double,3>& moment_nn = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM,2);
     
            for(unsigned int k=0; k<TDim; k++)
            {
                rData.U(i,k+1)   = moment[k];
                rData.Un(i,k+1)  = moment_n[k];
                rData.Unn(i,k+1) = moment_nn[k];
                //rData.f_ext(i,k)   = body_force[k];
            }
            rData.U(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            rData.Un(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY,1);
            rData.Unn(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY,2);
            
            rData.U(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY);
            rData.Un(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY,1);
            rData.Unn(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY,2);
            
            //rData.r(i) = this->GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
         }

    }

    //~ template< unsigned int TDim, unsigned int TNumNodes=TDim+1>
    double ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim>& DN_DX)
    {
        double h=0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

    // 3D tetrahedra shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,4,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    // 2D triangle shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }

    // 3D tetrahedra shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.25; Ncontainer(0,1) = 0.25; Ncontainer(0,2) = 0.25; Ncontainer(0,3) = 0.25;
    }

    // 2D triangle shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,3>& Ncontainer)
    {
        Ncontainer(0,0) = 1.0/3.0; Ncontainer(0,1) = 1.0/3.0; Ncontainer(0,2) = 1.0/3.0;
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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}

    ///@name Private Operations
    ///@{


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

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

} // namespace Kratos.

#endif // KRATOS_STOKES_ELEMENT_SYMBOLIC_3D_INCLUDED  defined
