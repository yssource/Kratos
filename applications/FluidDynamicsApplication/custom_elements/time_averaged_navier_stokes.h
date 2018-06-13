//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_TIME_AVERAGED_NAVIER_STOKES)
#define  KRATOS_TIME_AVERAGED_NAVIER_STOKES

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_elements/navier_stokes.h"
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

// TODO: UPDATE CLASS INFORMATION
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
class TimeAveragedNavierStokes : public NavierStokes<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef NavierStokes<TDim, TNumNodes>                     BaseType;
    typedef typename BaseType::MatrixType                   MatrixType;
    typedef typename BaseType::VectorType                   VectorType;
    typedef typename BaseType::IndexType                     IndexType;
    typedef typename BaseType::GeometryType               GeometryType;
    typedef typename BaseType::NodesArrayType           NodesArrayType;
    typedef typename BaseType::PropertiesType           PropertiesType;
    typedef typename BaseType::ElementDataStruct       ElementDataType;
    
    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(TimeAveragedNavierStokes);

    struct TimeAveragedElementDataStruct : public ElementDataType
    {
        BoundedMatrix<double, TNumNodes, TDim> vnnn;
        array_1d<double,TNumNodes> pnnn;

        unsigned int step;    // Current time-step
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TimeAveragedNavierStokes(
        IndexType NewId, 
        typename GeometryType::Pointer pGeometry)
        : NavierStokes<TDim, TNumNodes>(NewId, pGeometry){
    }

    TimeAveragedNavierStokes(
        IndexType NewId, 
        typename GeometryType::Pointer pGeometry, 
        typename PropertiesType::Pointer pProperties)
        : NavierStokes<TDim, TNumNodes>(NewId, pGeometry, pProperties){

    }

    /// Destructor.
    ~TimeAveragedNavierStokes(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& rThisNodes, 
        typename PropertiesType::Pointer pProperties) const override {

        KRATOS_TRY

        return Kratos::make_shared< TimeAveragedNavierStokes < TDim, TNumNodes > >(
            NewId, 
            this->GetGeometry().Create(rThisNodes), 
            pProperties);
            
        KRATOS_CATCH("");
    }

    Element::Pointer Create(
        IndexType NewId, 
        typename GeometryType::Pointer pGeom, 
        typename PropertiesType::Pointer pProperties) const override {

        KRATOS_TRY

        return Kratos::make_shared< TimeAveragedNavierStokes < TDim, TNumNodes > >(
            NewId, 
            pGeom, 
            pProperties);

        KRATOS_CATCH("");
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override {
            
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        TimeAveragedElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Gauss point position
        BoundedMatrix<double,TNumNodes, TNumNodes> N_container;
        this->GetShapeFunctionsOnGauss(N_container);

        for(unsigned int igauss = 0; igauss<N_container.size2(); igauss++)
        {
            noalias(data.N) = row(N_container, igauss);

            this->ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            // here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix  *= data.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in time-averaged Navier-Stokes element")
    }


    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override {

        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        TimeAveragedElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Gauss point position
        BoundedMatrix<double,TNumNodes, TNumNodes> N_container;
        this->GetShapeFunctionsOnGauss(N_container);

        // Loop on gauss point
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<N_container.size2(); igauss++)
        {
            noalias(data.N) = row(N_container, igauss);

            this->ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rRightHandSideVector) += rhs_local;
        }

        // rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

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
    std::string Info() const override {
        return "TimeAveragedNavierStokes #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info() << this->Id();
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

    ///@}
    ///@name Protected Operators
    ///@{

    TimeAveragedNavierStokes() : NavierStokes<TDim, TNumNodes>() {}

    ///@}
    ///@name Protected Operations
    ///@{

    void ComputeGaussPointLHSContribution(
        BoundedMatrix<double, TNumNodes *(TDim + 1),TNumNodes *(TDim + 1)> &lhs,
        const TimeAveragedElementDataStruct &data);

    void ComputeGaussPointRHSContribution(
        array_1d<double, TNumNodes *(TDim + 1)> &rhs,
        const TimeAveragedElementDataStruct &data);

    // Auxiliar function to fill the element data structure
    void FillElementData(
        TimeAveragedElementDataStruct& rData, 
        const ProcessInfo& rCurrentProcessInfo){

        // Base element fill element data call
        BaseType::FillElementData(rData, rCurrentProcessInfo);

        // Specific time-averaged element data fill
        rData.step = rCurrentProcessInfo[STEP];

        for (unsigned int i = 0; i < TNumNodes; ++i){
            rData.pnnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,3);
            const array_1d<double,3>& r_vavg_nnn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,3);
            for (unsigned int j = 0; j < TDim; ++j){
                rData.vnnn(i,j) = r_vavg_nnn[j];
            }
        }
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
