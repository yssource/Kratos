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
    // KRATOS_WATCH(mThisParameters["constraints_process_list"].size())

    // KRATOS_WATCH(mThisParameters["constraints_process_list"][4]["Parameters"]["interval"][1].GetDouble())
    double time = mrModelPart.GetProcessInfo()[TIME];
    // KRATOS_WATCH(time)
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::Execute()
{
    auto& process_info = mrModelPart.GetProcessInfo();
    bool cycle_found = false;
    std::vector<int> cycle_identificator;
    std::vector<int> cycles_after_advance_strategy;
    std::vector<double> damage;
    bool damage_indicator = false;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;
    bool cycles_from_last_advance = false;

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(DAMAGE, damage, process_info);
        r_elem.GetValueOnIntegrationPoints(CYCLES_AFTER_ADVANCE_STRATEGY, cycles_after_advance_strategy, process_info);
        KRATOS_WATCH(cycles_after_advance_strategy)
		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
                if (damage[i] > 0.0){
                    // KRATOS_WATCH(r_elem.Id())
                    // KRATOS_WATCH(i)
                    damage_indicator = true;
                }
                if (cycles_after_advance_strategy[i] > 1){
                    // KRATOS_WATCH(r_elem.Id())
                    // KRATOS_WATCH(i)
                    cycles_from_last_advance = true;
                }
        }
    }

    this->CyclePeriodPerIntegrationPoint(cycle_found);  //This method detects if a cycle has finished somewhere in the model and 
                                                        //computes the time period of the cycle that has just finished.
    
    if (cycle_found) {  //If a cycle has finished then it is possible to apply the advancing strategy
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy);  //Check if the conditions are optimal to apply the advancing strategy in 
                                                                        //terms of max stress and reversion factor variation.
        // KRATOS_WATCH(advancing_strategy);
        // KRATOS_WATCH(damage_indicator);       
        double increment = 0.0;
        // if (advancing_strategy & !damage_indicator & cycles_from_last_advance) {
        if (advancing_strategy & !damage_indicator) {
			this->TimeIncrement(increment);
			// KRATOS_WATCH(increment);
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
                    KRATOS_WATCH(r_elem.Id())
                    KRATOS_WATCH(i)
                }
        }

        r_elem.SetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.SetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        // KRATOS_WATCH(previous_cycle_time)
        // KRATOS_WATCH(period)

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
    bool fatigue_in_course = false;
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
            if (max_stress[i] > s_th[i]) {
		        std::cout << "WHAT \n";
                KRATOS_WATCH(r_elem.Id())
                fatigue_in_course = true;
                acumulated_max_stress_rel_error += max_stress_rel_error[i];
                acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
                // KRATOS_WATCH(i);
            }
        }
        // KRATOS_WATCH(acumulated_max_stress_rel_error)
        // KRATOS_WATCH(acumulated_rev_factor_rel_error)
    }
    if (acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4 && fatigue_in_course) {
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
    const double period_json = mThisParameters["fatigue"]["period"].GetDouble();
    bool current_constraints_process_list_detected = false;

    double user_avancing_time = mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble();
    double user_avancing_cycles = mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble();
    std::vector<std::string> constraints_list = mThisParameters["fatigue"]["constraints_process_list"].GetStringArray();
    
    // KRATOS_WATCH(constraints_list.size())
    // KRATOS_WATCH(constraints_list[0])
    // KRATOS_WATCH(mThisParameters["constraints_process_list"][0]["Parameters"]["model_part_name"].GetString())
    // KRATOS_WATCH(mThisParameters["constraints_process_list"][0]["Parameters"]["interval"][1].GetDouble())

    double model_part_final_time = mThisParameters["problem_data"]["end_time"].GetDouble();
    for (unsigned int i = 0; i < constraints_list.size(); i++){
        for (unsigned int j = 0; j < mThisParameters["constraints_process_list"].size(); j++){
            std::string model_part_name = mThisParameters["constraints_process_list"][j]["Parameters"]["model_part_name"].GetString();
            double model_part_end_time = mThisParameters["constraints_process_list"][j]["Parameters"]["interval"][1].GetDouble();
            KRATOS_WATCH(model_part_name)
            KRATOS_WATCH(model_part_end_time)
            KRATOS_WATCH(constraints_list[i])
            if (constraints_list[i] == model_part_name && time <= model_part_end_time && !current_constraints_process_list_detected) {
                model_part_final_time = model_part_end_time;
                current_constraints_process_list_detected = true;
                KRATOS_WATCH("no more")
            }
        }
    }
    double model_part_time_increment = model_part_final_time - time;
    // KRATOS_WATCH(model_part_final_time)
    // KRATOS_WATCH(model_part_time_increment)
    std::cout.precision(15); 
    std::cout << model_part_time_increment << "\n";
    // KRATOS_WATCH(time)


    // KRATOS_WATCH(mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble())
    // KRATOS_WATCH(mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble())
    // KRATOS_WATCH(mThisParameters["constraints_process_list"][4]["Parameters"]["interval"][1].GetDouble())
    // KRATOS_WATCH(user_avancing_time)
    // KRATOS_WATCH(user_avancing_cycles)
    // KRATOS_WATCH(mThisParameters["fatigue"]["constraints_process_list"].GetStringArray())

    min_time_to_failure = std::min(user_avancing_time, model_part_time_increment);
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info); 
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        r_elem.GetValueOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS, max_stress, process_info);
        KRATOS_WATCH(cycles_to_failure_element)
		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            if (max_stress[i] > s_th[i]) {
                // double Nf_time_to_failure = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
                double Nf_time_to_failure = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period_json;
                // double user_avancing_cycles_to_time = user_avancing_cycles * period[i];
                double user_avancing_cycles_to_time = user_avancing_cycles * period_json;
                if (Nf_time_to_failure < min_time_to_failure){
                    min_time_to_failure = Nf_time_to_failure;
                    // KRATOS_WATCH(min_time_to_failure)
                }
                if (user_avancing_cycles_to_time < min_time_to_failure){
                    min_time_to_failure = user_avancing_cycles_to_time;
                    // KRATOS_WATCH(min_time_to_failure)
                }
            }
        }
    }
	rIncrement = min_time_to_failure;

    KRATOS_WATCH(min_time_to_failure)
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::TimeAndCyclesUpdate(double Increment) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    std::vector<int>  local_number_of_cycles;
    std::vector<double> period;
    const double period_json = mThisParameters["fatigue"]["period"].GetDouble();

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, r_process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            //unsigned int local_cycles_increment;
            // if (period[i] == 0.0) {
            //     local_cycles_increment = 0;
            // } else {
            //     local_cycles_increment = std::trunc(Increment / period[i]);
            // }
            unsigned int local_cycles_increment = std::trunc(Increment / period_json);
            local_number_of_cycles[i] += local_cycles_increment;
        }
        KRATOS_WATCH(local_number_of_cycles)
        r_elem.SetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
    }
    r_process_info[TIME_INCREMENT] = std::trunc(Increment / period_json) * period_json;
    KRATOS_WATCH(Increment)

    // time += Increment;
    // r_process_info[TIME] = time;
}

} // namespace Kratos
