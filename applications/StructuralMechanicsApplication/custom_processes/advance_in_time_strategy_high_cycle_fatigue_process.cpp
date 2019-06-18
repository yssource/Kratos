//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#include "custom_processes/advance_in_time_strategy_high_cycle_fatigue_process.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

// Constructor
	AdvanceInTimeStrategyHighCycleFatigueProcess::AdvanceInTimeStrategyHighCycleFatigueProcess(
		ModelPart &rModelPart,
		Parameters ThisParameters) : mrModelPart(rModelPart), mThisParameters(ThisParameters)
{
    KRATOS_WATCH(mThisParameters["constraints_process_list"].size())
    KRATOS_WATCH(mThisParameters["constraints_process_list"][4]["Parameters"]["interval"][1].GetDouble())
    double time = mrModelPart.GetProcessInfo()[TIME];
    KRATOS_WATCH(time)
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::Execute()
{
    auto& process_info = mrModelPart.GetProcessInfo();
    bool cycle_found = false;
    std::vector<int> cycle_identificator;
    std::vector<double> damage;
    bool damage_indicator = false;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(DAMAGE, damage, process_info);
		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
                if (damage[i] > 0.0){
                    damage_indicator = true;
                }
        }
    }

    this->CyclePeriodPerIntegrationPoint(cycle_found);
    
    if (cycle_found) {
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy);
        KRATOS_WATCH(advancing_strategy);
        KRATOS_WATCH(damage_indicator);        

        double increment = 0.0;
        if (advancing_strategy & !damage_indicator) {
			this->TimeIncrement(increment);
			KRATOS_WATCH(increment);
            if(increment > 0){
                this->TimeAndCyclesUpdate(increment);
                process_info[ADVANCE_STRATEGY_APPLIED] = true;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::CyclePeriodPerIntegrationPoint(bool& rCycleFound) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<bool> cycle_identificator;
    std::vector<double> previous_cycle_time;
    std::vector<double> period;
    double time = process_info[TIME];
    
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLE_INDICATOR, cycle_identificator, process_info);
        r_elem.GetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        KRATOS_WATCH(previous_cycle_time)
        KRATOS_WATCH(period)


		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
                if (cycle_identificator[i]){
                    period[i] = time - previous_cycle_time[i];
                    previous_cycle_time[i] = time;
                    rCycleFound = true;
                }
        }

        r_elem.SetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.SetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        KRATOS_WATCH(previous_cycle_time)
        KRATOS_WATCH(period)

    }

}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::StableConditionForAdvancingStrategy(bool& rAdvancingStrategy) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double> max_stress_rel_error;
    std::vector<double> rev_factor_rel_error;
    
    double acumulated_max_stress_rel_error;
    double acumulated_rev_factor_rel_error;
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.Id();
        acumulated_max_stress_rel_error = 0.0;
        acumulated_rev_factor_rel_error = 0.0;
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);

        // double elemento = r_elem;
        // KRATOS_WATCH(r_elem);
        // KRATOS_WATCH(elemento);
        for (unsigned int i=0; i < max_stress_rel_error.size(); i++) {
            acumulated_max_stress_rel_error += max_stress_rel_error[i];
            acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
        }
        if (acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4) {
            rAdvancingStrategy = true;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::TimeIncrement(double& rIncrement) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double>  cycles_to_failure_element;
    std::vector<int>  local_number_of_cycles;
    std::vector<double> period;
    double min_time_to_failure = 0.0;

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info); 
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            double time_to_failure = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
            if (period[i] > 0 || min_time_to_failure == 0.0 || time_to_failure < min_time_to_failure){
                min_time_to_failure = time_to_failure;
                KRATOS_WATCH(min_time_to_failure)
            }
        }
    }
	rIncrement = min_time_to_failure;
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::TimeAndCyclesUpdate(double Increment) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    std::vector<int>  local_number_of_cycles;
    std::vector<double> period;

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, r_process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            unsigned int local_cycles_increment = std::trunc(Increment / period[i]);
            local_number_of_cycles[i] += local_cycles_increment;
        }
        r_elem.SetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
    }
    // double time = r_process_info[TIME];
    // time += Increment;
    // r_process_info[TIME] = time;
}

} // namespace Kratos
