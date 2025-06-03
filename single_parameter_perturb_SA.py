"""
Sensitivity Analysis

sensitivity_N

Author: Halima sadia
Date: 2025-01-17
"""


import numpy as np
from equations_SA import ODEsystem_SA
from scipy.integrate import solve_ivp 
import parameters as param
from InitialConditions import InitialConditions
import matplotlib.pyplot as plt
AgeStart = 30
AgeEnd = 80
t_eval = np.linspace(365 * AgeStart, 365 * AgeEnd, 1000)
# run the model to solve the ODES 
def run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_options):
    y0 = InitialConditions(p, AgeStart)
    try:
        sol = solve_ivp(ODEsystem_SA, [365 * AgeStart, 365 * AgeEnd], y0, method="BDF", args=[p], t_eval=t_eval, **solver_options)
        if not sol.success:
            print(f"Integration failed: {sol.message}")
            return None
        return sol
    except Exception as e:
        print(f"An error occurred during integration: {str(e)}")
        return None


cases =[
    (0, 0, "Women (APOE-)"),
    (0, 1, "Women (APOE+)"),
    (1, 0, "Men (APOE-)"),
    (1, 1, "Men (APOE+)")
]  
# originl solution
def perturbation_parameter_plot(case, perturbation_ratio, solver_options):
    sex, apoe4_status, case_name=case
    p = param.Parameters(sex, apoe4_status) 
    parameter_name = 'd_Fi'
    original_value = getattr(p, parameter_name)
    print(original_value)


    # Create a figure for Neuron count (N)
    fig, axs = plt.subplots(1, 1, figsize=(14, 12), dpi=250)
    fig.suptitle(f'Parameter Perturbations for {case_name} - Neuron Count',  fontsize=24, y=0.04) 

    # loop through the perturbation_ratio
    # perturb solution 
    for perturbation_ratio in [1] + perturbation_ratio:
        p = param.Parameters(sex, apoe4_status)
        setattr(p, parameter_name, perturbation_ratio * original_value )
        sol = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_options)

        if sol is not None:
            
         axs.plot(sol.t / 365.25, sol.y[8], label=f'{perturbation_ratio*100-100:+.0f}%')
         axs.set_title(f'{parameter_name}, - Neuron_count(N)', fontsize=24)
         axs.set_xlabel('Time(years)', fontsize=26)
         axs.set_ylabel('Neuron_count(N)', fontsize=26)
         axs.legend(fontsize=18)
    
    #plt.figtext(0.5, 0.01, f'Parameter Perturbations for {case_name} - Neuron Count', ha='center', fontsize=14)

    plt.tight_layout()
    plt.subplots_adjust(top=0.96, bottom=0.15, hspace=0.4)   
    plt.savefig(f'parameter_perturbations_N_{case_name.replace(" ", "_")}.png', dpi=250)
    plt.show()
perturbation_ratio = [1.05, 1.10, 0.95]  # 5%, 10%, and -5% perturbation     
solver_options = {'atol': 1e-10, 'rtol': 1e-10}
    # Run the perturbation analysis and plot the results for each case
for case in cases:
    print(f"Processing case: {case[2]}")
    perturbation_parameter_plot(case, perturbation_ratio, solver_options)

print("Analysis completed for all cases.")


     