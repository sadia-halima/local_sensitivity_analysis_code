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

# Define the age range for the simulation
AgeStart = 30
AgeEnd = 80
t_eval = np.linspace(365 * AgeStart, 365 * AgeEnd, 1000)

def run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_option):
    y0 = InitialConditions(p, AgeStart)
    try:
        sol = solve_ivp(ODEsystem_SA, [365 * AgeStart, 365 * AgeEnd], y0, method="BDF", args=[p], t_eval=t_eval, **solver_option)
        if not sol.success:
            print(f"Integration failed: {sol.message}")
            return None
        return sol
    except Exception as e:
        print(f"An error occurred during integration: {str(e)}")
        return None

def perturb_parameters_and_plot(case, perturbation_factors, solver_options):
    sex, apoe4_status, case_name = case
    parameters_to_analyze = ['d_Fi', 'lambda_Gtau']

    for parameter_name in parameters_to_analyze:
        p = param.Parameters(sex, apoe4_status)
        if hasattr(p, parameter_name):
            original_value = getattr(p, parameter_name)

            # Define biomarker labels and their corresponding indices
            biomarkers = [
    ('Amyloid-beta', 3, 'μM'),         # micromolar
    ('Neuron Count', 8, 'neurons/mL'),  # neurons per milliliter
    ('Tau', 6, 'μM')                    # micromolar
]

            # Loop through each biomarker and create a separate figure
            for biomarker, index, unit in biomarkers:
                plt.figure(figsize=(10, 6), dpi=250)
                plt.title(f'{parameter_name} Perturbations for {case_name} - {biomarker}', fontsize=18)

                # Loop through perturbation factors
                for perturbation_factor in [1] + perturbation_factors:
                    p = param.Parameters(sex, apoe4_status)
                    setattr(p, parameter_name, original_value * perturbation_factor)
                    sol = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_options)

                    if sol is not None:
                        plt.plot(sol.t / 365.25, sol.y[index], label=f'{perturbation_factor*100-100:+.0f}%')

                plt.xlabel('Age (years)', fontsize=16)
                plt.ylabel(f'{biomarker} [$\mathrm{{{unit}}}$]', fontsize=16)
                plt.legend(fontsize=12)
                plt.grid(True)
                plt.savefig(f'{parameter_name}_perturbations_{case_name.replace(" ", "_")}_{biomarker.replace(" ", "_")}.jpg', dpi=300)
                plt.show()  # Show each figure separately

        else:
            print(f"Parameter '{parameter_name}' not found in the model parameters.")

# Set solver tolerances
solver_options = {'atol': 1e-22, 'rtol': 1e-5}

# Define cases and perturbation factors
cases = [
    (0, 0, "Women (APOE-)"),
    #(0, 1, "Women (APOE+)"),
    #(1, 0, "Men (APOE-)"),
    #(1, 1, "Men (APOE+)")
]
perturbation_factors = [1.05, 1.10, 0.95]  # 5%, 10%, and -5% perturbation

# Run the perturbation analysis and plot the results for each case
for case in cases:
    print(f"Processing case: {case[2]}")
    perturb_parameters_and_plot(case, perturbation_factors, solver_options)

print("Analysis completed for all cases.")
