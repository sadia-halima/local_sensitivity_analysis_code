# local_sensitivity_analysis_code
Code for Sensitivity Analysis of a mathematical model of AD in the submitted manuscript
This repository contains the code for sensitivity analysis of a mathematical model of Alzheimer's disease progression. The model incorporates the dynamics of amyloid-beta aggregation, tau protein phosphorylation, neuron degeneration, and immune response, including microglia and macrophages.

The purpose of this sensitivity analysis is to identify the most influential parameters driving changes in key disease biomarkers such as:
- Amyloid-beta (Aβ)
- Tau protein (τ)
- Living neurons (N)

## 📁 Repository Structure

Alzheimers_Model/
│
├── src/
│   ├── equations_SA.py                          # System of 19 ODEs representing biological processes

│   ├── InitialConditions.py                     # Computes initial values for all variables
│   ├── parameters.py                            # Parameter definitions (must be added)
│
├── analysis/
│   ├── mean_relative_change.py                  # Sensitivity using mean relative % change
│   ├── One-at-a-time_SA.py                      # One-at-a-time parameter perturbation
│   ├── groupwiseSynergy.py                      # Group-wise sensitivity or synergy exploration
│   ├── specific_parameter_perturbation_OAT.py   # Specific OAT perturbation for selected parameters
│   ├── single_parameter_perturb.py              # Single parameter test analysis
│
├── results/
│   ├── mean_relative_sensitivity_AB.png         # Aβ sensitivity results
│   ├── mean_relative_sensitivity_N.png          # Neuron loss sensitivity results
 mean_relative_sensitivity_tau.png        # Tau sensitivity results

requirements.txt                             # Python package dependencies
 README.md                                    # This documentation


## 🚀 How to Run the Code

### 1. Clone the Repository
https://github.com/sadia-halima/local_sensitivity_analysis_code.git
```bash

### 2. pip install -r requirements.txt
numpy
pandas
scipy
matplotlib

