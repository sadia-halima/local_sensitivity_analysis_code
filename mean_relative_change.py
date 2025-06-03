"""This file contains the equations of the model.
It defines a function that generates the ODE system of the model.

:Author: Halima Sadia
:Creation date: oct 5, 2024
:Last modified:


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import parameters as param
from InitialConditions import InitialConditions
from equations_SA import ODEsystem_SA

# Age parameters
AgeStart, AgeEnd = 30, 80
t_eval = np.linspace(365 * AgeStart, 365 * AgeEnd, 1000)

# Solver tolerances
solver_options = {'atol': 1e-10, 'rtol': 1e-10} 


# Run model function
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

def relative_change(Sex, APOE4_status, perturbation_factor, var_int, solver_option):
    p = param.Parameters(Sex, APOE4_status)
    AgeStart, AgeEnd = 30, 80
    t_eval = np.linspace(365 * AgeStart, 365 * AgeEnd, 1000)

    sol_original = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_option)
    if sol_original is None:
        print("Original model integration failed")
        return pd.DataFrame()

    age_measure = -1
    original_output = sol_original.y[var_int][age_measure]
    relative_changes = []

    # Iterate over parameters, excluding "S" (Sex) and "AP" (APOE4 status)
    for attr in p.__dict__:
        if attr in ["S", "AP"]:  # Skip these parameters
            continue

        p = param.Parameters(Sex, APOE4_status)
        original_value = getattr(p, attr)

        # Perturb positively
        p_increased = original_value * perturbation_factor
        setattr(p, attr, p_increased)
        sol_increased = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_option)
        if sol_increased is None:
            print(f"Skipping parameter {attr} due to integration failure (positive perturbation)")
            continue
        relative_change_increase = np.abs((sol_increased.y[var_int][age_measure] - original_output) / original_output)

        # Perturb negatively
        p_decreased = original_value / perturbation_factor
        setattr(p, attr, p_decreased)
        sol_decreased = run_model_SA(p, AgeStart, AgeEnd, t_eval, solver_option)
        if sol_decreased is None:
            print(f"Skipping parameter {attr} due to integration failure (negative perturbation)")
            continue
        relative_change_decrease = np.abs((sol_decreased.y[var_int][age_measure] - original_output) / original_output)

        # Compute the mean of both changes
        relative_change = 100 * np.mean([relative_change_increase, relative_change_decrease])

        relative_changes.append({
            "Parameter": attr,
            "Relative Change": relative_change
        })

        # Restore original value
        setattr(p, attr, original_value)

    df_relative_changes = pd.DataFrame(relative_changes)
    return df_relative_changes.sort_values(by="Relative Change", ascending=False)


# Cases for analysis
cases = [
    (0, 0, "Women (APOE-)"),
    (0, 1, "Women (APOE+)"),
    (1, 0, "Men (APOE-)"),
    (1, 1, "Men (APOE+)")
]

# Store results for each case
all_results = {'AB': {}, 'N': {}, 'tau': {}}

for sex, apoe_status, title in cases:
    # Calculate relative changes for N, AB, and Tau
    all_results['N'][title] = relative_change(sex, apoe_status, 1.1, 8, solver_options)
    all_results['AB'][title] = relative_change(sex, apoe_status, 1.1, 3, solver_options)
    all_results['tau'][title] = relative_change(sex, apoe_status, 1.1, 6, solver_options)



def plot_relative_change(biomarker, filename):
    data = all_results[biomarker]
    
    # Get union of all parameters across all cases
    all_params = set()
    for df in data.values():
        all_params.update(df['Parameter'])
    all_params = sorted(list(all_params))
    
    # Create a DataFrame with all parameters and cases
    df = pd.DataFrame(index=all_params, columns=[case[2] for case in cases])
    
    # Fill the DataFrame with mean relative changes
    for case in cases:
        title = case[2]
        df[title] = data[title].set_index('Parameter')['Relative Change']

    # Handle missing parameters (e.g., N_0) for specific biomarkers
    if biomarker == 'N':
        df = df[df.index != 'N_0']

    # Apply threshold to filter out less significant parameters
    highest_relative_change = df.max().max()
    threshold = 0.2 * highest_relative_change
    df = df[df.max(axis=1) >= threshold]
    
    # Sort parameters by total sensitivity
    df = df.fillna(0)
    df['Total'] = df.mean(axis=1)  # Calculate mean relative change across all cases
    df = df.sort_values('Total', ascending=False)
    df = df.drop('Total', axis=1)  # Remove the auxiliary 'Total' column

    # Debugging: Check DataFrame
    print(f"DataFrame for {biomarker}:\n", df)

    # Skip empty DataFrame
    if df.empty:
        print(f"No data available for biomarker {biomarker}, skipping plot.")
        return

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(22, 14))
    bar_width = 0.2
    index = np.arange(len(df))
    
    for i, case in enumerate(cases):
        ax.bar(index + i * bar_width, df[case[2]], bar_width, label=case[2])
    
    # Set labels and title
    ax.set_xticks(index + bar_width * len(cases) / 2)
    ax.set_xticklabels(df.index, rotation=45, ha='right', fontsize=26)
    ax.set_xlabel('Parameters', fontsize=26)  # X-label
    ax.set_ylabel(f'Relative Change of {biomarker} (%)',  fontsize=26)
    ax.set_title(f'Relative Change of {biomarker} in response to a 10% change in parameter', fontsize=26)
    ax.legend(title='Cases', loc='upper right', fontsize=22, bbox_to_anchor=(1, 1))
    
    # Adjust layout
    plt.tight_layout(pad=9.5)
    plt.subplots_adjust(bottom=0.2, left=0.1, right=0.85, top=0.9)
    fig.autofmt_xdate(rotation=45)
    
    # Debugging: Check layout
    plt.gcf().canvas.draw()  # Force redraw for debugging
    
    # Save and display plot
    plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()


# Plot and save the figures for AB, N, and Tau
plot_relative_change('AB', 'mean_relative_sensitivity_AB.png')
plot_relative_change('N', 'mean_relative_sensitivity_N.png')
plot_relative_change('tau', 'mean_relative_sensitivity_tau.png')

print("All figures with mean relative changes have been saved.")

