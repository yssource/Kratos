//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Author:          Jordi Cotela
//


#ifndef KRATOS_TURBULENCE_STATISTICS_PROCESS_H
#define KRATOS_TURBULENCE_STATISTICS_PROCESS_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/turbulence_statistics_container.h"


namespace Kratos
{

///@name Kratos Classes
///@{

class TurbulenceStatisticsProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceStatisticsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TurbulenceStatisticsProcess(ModelPart::Pointer pModelPart,
                                double StartTime,
                                bool ResetContainers):
        mpModelPart(pModelPart),
        mStartTime(StartTime),
        mResetContainers(ResetContainers)
    {}

    /// Destructor.
    virtual ~TurbulenceStatisticsProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    virtual void Execute() override
    {
        KRATOS_TRY;
        KRATOS_THROW_ERROR(std::runtime_error,"Do not call TurbulenceStatisticsProcess::Execute. It does nothing","");
        KRATOS_CATCH("");
    }


    /**
     * @brief ExecuteBeforeSolutionLoop Initialize containers.
     */
    virtual void ExecuteBeforeSolutionLoop() override
    {
        // Skip this step if we are reading from a restart file and we already recorded results before
        if ( mResetContainers )
        {
            mpModelPart->GetProcessInfo()[RECORDED_STEPS] = 0;

//            for ( ModelPart::NodeIterator i = mpModelPart->GetCommunicator().LocalMesh().NodesBegin();
//                  i != mpModelPart->GetCommunicator().LocalMesh().NodesEnd(); i++)
//            {
//                i->GetValue(MEAN_VELOCITY) = array_1d<double,3>(3,0.0);
//                i->GetValue(MEAN_PRESSURE) = 0.0;
//                i->GetValue(VELOCITY_COVARIANCES).resize(3,3);
//                noalias( i->GetValue(VELOCITY_COVARIANCES) ) = ZeroMatrix(3,3);
//            }
            //for ( ModelPart::ElementIterator i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin();
            //      i != mpModelPart->GetCommunicator().LocalMesh().ElementsEnd(); i++)
            for (int j = 0; j < static_cast<int>(mpModelPart->GetCommunicator().LocalMesh().Elements().size()); j++)
            {
                auto i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin() + j;
                unsigned int NumGauss = i->GetGeometry().IntegrationPointsNumber(i->GetIntegrationMethod());
                i->GetValue(TURBULENCE_STATISTICS) = TurbulenceStatisticsContainer::Pointer( new TurbulenceStatisticsContainer(NumGauss) );
            }
        }
    }


    virtual void ExecuteInitializeSolutionStep() override
    {

    }

    /**
     * @brief ExecuteFinalizeSolutionStep Add solution step results to measured averages.
     */
    virtual void ExecuteFinalizeSolutionStep() override
    {
        if ( mpModelPart->GetProcessInfo()[TIME] >= mStartTime)
        {
            mpModelPart->GetProcessInfo()[RECORDED_STEPS]++;

//            for ( ModelPart::NodeIterator i = mpModelPart->GetCommunicator().LocalMesh().NodesBegin();
//                  i != mpModelPart->GetCommunicator().LocalMesh().NodesEnd(); i++)
//            {
//                const array_1d<double,3>& rU = i->FastGetSolutionStepValue(VELOCITY);

//                i->GetValue(MEAN_VELOCITY) += rU;
//                i->GetValue(VELOCITY_COVARIANCES) += boost::numeric::ublas::outer_prod(rU,rU);

//                i->GetValue(MEAN_PRESSURE) += i->FastGetSolutionStepValue(PRESSURE);
//            }

            std::vector<double> out;
            ProcessInfo& rProcessInfo = mpModelPart->GetProcessInfo();
            //for ( ModelPart::ElementIterator i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin();
            //      i != mpModelPart->GetCommunicator().LocalMesh().ElementsEnd(); i++)
            for (int j = 0; j < static_cast<int>(mpModelPart->GetCommunicator().LocalMesh().Elements().size()); j++)
            {
                auto i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin() + j;
                i->GetValueOnIntegrationPoints(MEAN_KINETIC_ENERGY,out,rProcessInfo);
            }
        }
    }

//    /**
//     * @brief ExecuteBeforeOutputStep Compute partial results and set them for printig through GiDIO.
//     * @note This function requires having historical nodal data containers for nodal statistics!
//     */
//    virtual void ExecuteBeforeOutputStep()
//    {
//        if ( mpModelPart->GetProcessInfo()[TIME] >= mStartTime )
//        {
//            double NumSteps = double(mpModelPart->GetProcessInfo()[RECORDED_STEPS]);

//            for ( ModelPart::NodeIterator i = mpModelPart->GetCommunicator().LocalMesh().NodesBegin();
//                  i != mpModelPart->GetCommunicator().LocalMesh().NodesEnd(); i++)
//            {
//                i->FastGetSolutionStepValue(MEAN_VELOCITY) = i->GetValue(MEAN_VELOCITY) / NumSteps;
//                i->FastGetSolutionStepValue(MEAN_PRESSURE) = i->GetValue(MEAN_PRESSURE) / NumSteps;

//                // We use the relationship var(u) = mean(u^2) - mean(u)^2
//                i->FastGetSolutionStepValue(VELOCITY_COVARIANCES).resize(3,3,false);
//                Matrix& CVuj = i->FastGetSolutionStepValue(VELOCITY_COVARIANCES);

//                const Matrix& MeanUSquared = i->GetValue(VELOCITY_COVARIANCES);
//                noalias(CVuj) = MeanUSquared / NumSteps;
//                const array_1d<double,3>& MeanU = i->FastGetSolutionStepValue(MEAN_VELOCITY);

//                for (unsigned int m = 0; m < 3; m++)
//                    for (unsigned int n = 0; n < 3; n++)
//                        CVuj(m,n) -= MeanU[m]*MeanU[n];
//            }

//            TurbulenceStatisticsContainer::DumpData(*mpModelPart);
//        }
//    }

    /**
     * @brief ExecuteBeforeOutputStep Call this to generate a dump of several instantaneous values of interest
     */
    virtual void ExecuteBeforeOutputStep() override
    {
        std::vector<double> out;
        //for ( ModelPart::ElementIterator i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin();
        //      i != mpModelPart->GetCommunicator().LocalMesh().ElementsEnd(); i++)
        //{
        for (int j = 0; j < static_cast<int>(mpModelPart->GetCommunicator().LocalMesh().Elements().size()); j++)
        {
            auto i = mpModelPart->GetCommunicator().LocalMesh().ElementsBegin() + j;
            i->GetValueOnIntegrationPoints(TAU,out,mpModelPart->GetProcessInfo());
        }

        TurbulenceStatisticsContainer::DumpInstantData(*mpModelPart);
    }


    /**
     * @brief ExecuteFinalize Finalize measured quantities.
     */
    virtual void ExecuteFinalize() override
    {
        if ( mpModelPart->GetProcessInfo()[TIME] >= mStartTime)
        {
//            double NumSteps = double(mpModelPart->GetProcessInfo()[RECORDED_STEPS]);

//            for ( ModelPart::NodeIterator i = mpModelPart->GetCommunicator().LocalMesh().NodesBegin();
//                  i != mpModelPart->GetCommunicator().LocalMesh().NodesEnd(); i++)
//            {
//                // We use the relationship var(u) = mean(u^2) - mean(u)^2
//                Matrix& CVuj = i->GetValue(VELOCITY_COVARIANCES);
//                array_1d<double,3>& MeanU = i->GetValue(MEAN_VELOCITY);
//                for (unsigned int m = 0; m < 3; m++)
//                    for (unsigned int n = 0; n < 3; n++)
//                        CVuj(m,n) -= MeanU[m]*MeanU[n];

//                CVuj /= NumSteps;
//                MeanU /= NumSteps;
//                i->GetValue(MEAN_PRESSURE) /= NumSteps;
//            }

            TurbulenceStatisticsContainer::DumpData(*mpModelPart);
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
    virtual std::string Info() const override
    {
        return "TurbulenceStatisticsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TurbulenceStatisticsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{

    ModelPart::Pointer mpModelPart;

    const double mStartTime;

    const bool mResetContainers;

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TurbulenceStatisticsProcess& operator=(TurbulenceStatisticsProcess const& rOther);

    /// Copy constructor.
    TurbulenceStatisticsProcess(TurbulenceStatisticsProcess const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TurbulenceStatisticsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TurbulenceStatisticsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_TURBULENCE_STATISTICS_PROCESS_H
