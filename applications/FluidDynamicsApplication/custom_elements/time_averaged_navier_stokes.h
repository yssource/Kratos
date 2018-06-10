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

    typedef NavierStokes<TDim, TNumNodes>               BaseType;
    typedef typename BaseType::IndexType               IndexType;
    typedef typename BaseType::GeometryType         GeometryType;
    typedef typename BaseType::NodesArrayType     NodesArrayType;
    typedef typename BaseType::PropertiesType     PropertiesType;
    
    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(TimeAveragedNavierStokes);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, TDim> y, yn, ynn, ynnn, vmesh, f;
        array_1d<double,TNumNodes> x, xn, xnn, xnnn, rho, mu;

        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;

        Matrix C;
        Vector stress;
        Vector strain;

        double bdf0;
        double bdf1;
        double bdf2;
        double c;             // Wave velocity (needed if artificial compressibility is considered)
        double h;             // Element size
        double volume;        // In 2D: element area. In 3D: element volume
        double dt;            // Time increment
        double dyn_tau;       // Dynamic tau considered in ASGS stabilization coefficients

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
        const ElementDataStruct &data);

    void ComputeGaussPointRHSContribution(
        array_1d<double, TNumNodes *(TDim + 1)> &rhs,
        const ElementDataStruct &data);

    // Auxiliar function to fill the element data structure
    void FillElementData(
        ElementDataStruct& rData, 
        const ProcessInfo& rCurrentProcessInfo){

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

        rData.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];  // Only, needed if the temporal dependent term is considered in the subscales
        rData.dt = rCurrentProcessInfo[DELTA_TIME];        // Only, needed if the temporal dependent term is considered in the subscales

        rData.c = rCurrentProcessInfo[SOUND_VELOCITY];     // Wave velocity

        rData.step = rCurrentProcessInfo[STEP];            // Current time step

        for (unsigned int i = 0; i < TNumNodes; i++){

            const array_1d<double,3>& r_vavg = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& r_vavg_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& r_vavg_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& r_vavg_nnn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,3);
            const array_1d<double,3>& r_vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double,3>& r_body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);

            for(unsigned int k=0; k<TDim; k++){
                rData.y(i,k) = r_vavg[k];
                rData.yn(i,k) = r_vavg_n[k];
                rData.ynn(i,k) = r_vavg_nn[k];
                rData.ynnn(i,k) = r_vavg_nnn[k];
                rData.vmesh(i,k) = r_vel_mesh[k];
                rData.f(i,k) = r_body_force[k];
            }

            rData.x[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            rData.xn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,1);
            rData.xnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,2);
            rData.xnnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,3);
            rData.rho[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            rData.mu[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
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
