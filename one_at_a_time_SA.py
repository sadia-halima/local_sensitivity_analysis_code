"""
Sensitivity Analysis

sensitivity_testing_AB

Author: Halima Sadia
Date: feb 2025
"""

"""
General Explanations:

Axes: The horizontal axis (x-axis) tracks time (in days) of the simulation.
 The vertical axis (y-axis) shows the concentration of AB (N, Tau).

Lines: The graph displays two lines representing the AB (N, Tau) value over time. 
One line is for the original parameter value (original),
 and the other shows the AB (N, Tau) value when the parameter is slightly increased (perturbed) by 5%.

Parameter Sensitivity: Parameters that significantly alter the AB (N, Tau) value are considered more sensitive to changes.

Interpreting Effects:
Positive vs. Negative Effects: A parameter has a positive effect on AB (N, Tau) if 
increasing the parameter value leads to a higher AB (N, Tau) concentration.
 Conversely, a parameter has a negative effect if increasing its value results in a lower AB (N, Tau) concentration.
"""

import numpy as np
from scipy.integrate import solve_ivp
import parameters as param
from InitialConditions import InitialConditions
from equations_SA import ODEsystem_SA
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    p = param.Parameters(sex, apoe4_status)
    num_parameters = len(p.__dict__)
    num_perturbations = len(perturbation_factors) + 1  # Including original

    fig, axs = plt.subplots(num_parameters, 3, figsize=(20, num_parameters * 5), dpi=100)
    fig.suptitle(f'Parameter Perturbations for {case_name}', fontsize=16, y=0.964)  # Adjust y value

    for param_index, parameter_name in enumerate(p.__dict__.keys()):
        original_value = getattr(p, parameter_name)

        for perturb_index, perturbation_factor in enumerate([1] + perturbation_factors):
            p = param.Parameters(sex, apoe4_status)
            setattr(p, parameter_name, original_value * perturbation_factor)
            sol = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_options)

            if sol is not None:
                for i, (biomarker, index) in enumerate([('AB', 3), ('N', 8), ('tau', 6)]):
                    ax = axs[param_index, i]
                    ax.plot(sol.t / 365.25, sol.y[index], label=f'{perturbation_factor*100-100:+.0f}%')
                    ax.set_title(f'{parameter_name} - {biomarker}')
                    ax.set_xlabel('Time (years)')
                    ax.set_ylabel(f'{biomarker} Value')
                    ax.legend()

    plt.tight_layout()
    plt.subplots_adjust(top=0.96, hspace=0.4)   
    plt.savefig(f'parameter_perturbations_{case_name.replace(" ", "_")}.png', dpi=100)
    plt.show()


# Set solver tolerances
solver_options = {'atol': 1e-22, 'rtol': 1e-5}

# Define cases and perturbation factors
cases = [
    (0, 0, "Women (APOE-)"),
    (0, 1, "Women (APOE+)"),
    (1, 0, "Men (APOE-)"),
    (1, 1, "Men (APOE+)")
]
perturbation_factors = [1.05, 1.10, 0.95]  # 5%, 10%, and -5% perturbation

# Run the perturbation analysis and plot the results for each case
for case in cases:
    print(f"Processing case: {case[2]}")
    perturb_parameters_and_plot(case, perturbation_factors, solver_options)

print("Analysis completed for all cases.")
