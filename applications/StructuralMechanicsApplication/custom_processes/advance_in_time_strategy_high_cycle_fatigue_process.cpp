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

    this->CyclePeriodPerIntegrationPoint(cycle_found);  //This method detects if a cycle has finished somewhere in the model and 
                                                        //computes the time period of the cycle that has just finished.
    
    if (cycle_found) {  //If a cycle has finished then it is possible to apply the advancing strategy
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy);  //Check if the conditions are optimal to apply the advancing strategy in 
                                                                        //terms of max stress and reversion factor variation.
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
    std::vector<double> previous_cycle_time;    //time when the previous cycle finished. It is used to obtain the new period for the current cycle
    std::vector<double> period;
    double time = process_info[TIME];
    
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLE_INDICATOR, cycle_identificator, process_info);
        r_elem.GetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        // KRATOS_WATCH(previous_cycle_time)
        // KRATOS_WATCH(period)


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
    rAdvancingStrategy = false;
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double> max_stress_rel_error;
    std::vector<double> rev_factor_rel_error;
    std::vector<double> s_th;
    std::vector<double> max_stress;
    
    double acumulated_max_stress_rel_error = 0.0;
    double acumulated_rev_factor_rel_error = 0.0;
    for (auto& r_elem : mrModelPart.Elements()) {   //Este loop se hace por todos los elementos y PI independientemente que haya habido o no cambio 
                                                    //de ciclo, porque se debe garantizar condición de estabilidad en TODO el modelo (siempre que 
                                                    //Smax > Sth)
        
        r_elem.Id();
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

        KRATOS_WATCH(s_th);
        KRATOS_WATCH(max_stress);
        KRATOS_WATCH(rev_factor_rel_error)
        KRATOS_WATCH(max_stress_rel_error)

        

        // for (unsigned int i=0; i < max_stress_rel_error.size(); i++) {
        const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            if (max_stress[i] > s_th[1]) {
                acumulated_max_stress_rel_error += max_stress_rel_error[i];
                acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
                KRATOS_WATCH(i);
            }
        }
        KRATOS_WATCH(acumulated_max_stress_rel_error)
        KRATOS_WATCH(acumulated_rev_factor_rel_error)
    }
    if (acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4) {
		std::cout << "NNNNNNOOOOOO \n";
        rAdvancingStrategy = true;
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
    std::vector<double> s_th;
    std::vector<double> max_stress;
    double min_time_to_failure;
    double time = process_info[TIME];

    double user_avancing_time = mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble();
    double user_avancing_cycles = mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble();
    KRATOS_WATCH(mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble())
    KRATOS_WATCH(mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble())
    KRATOS_WATCH(mThisParameters["constraints_process_list"][4]["Parameters"]["interval"][1].GetDouble())
    KRATOS_WATCH(user_avancing_time)
    KRATOS_WATCH(user_avancing_cycles)

    min_time_to_failure = user_avancing_time;
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info); 
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        r_elem.GetValueOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            if (max_stress[i] > s_th[1]) {
                double Nf_time_to_failure = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
                double user_avancing_cycles_to_time = user_avancing_cycles * period[i];
                if (Nf_time_to_failure < min_time_to_failure){
                    min_time_to_failure = Nf_time_to_failure;
                    KRATOS_WATCH(min_time_to_failure)
                }
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
        KRATOS_WATCH(local_number_of_cycles)
        r_elem.SetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
    }
    // double time = r_process_info[TIME];
    // time += Increment;
    // r_process_info[TIME] = time;
}

} // namespace Kratos
