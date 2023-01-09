import sys
BMS_FOLDER = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS"
sys.path.append(BMS_FOLDER)
sys.path.append(BMS_FOLDER + "\\" + "Prior")
from utils.kinetic_utils import PropeneAmmoxidation
import numpy as np
import os

OPS_mono = {
    'sin': 1,
    'cos': 1,
    'tan': 1,
    'exp': 1,
    'log': 1,
    'sinh' : 1,
    'cosh' : 1,
    'tanh' : 1,
    'pow2' : 1,
    'pow3' : 1,
#    'abs'  : 1,
    'sqrt' : 1,
#    'fac' : 1,
    '-' : 1,
    '+' : 2,
    '*' : 2,
    '/' : 2,
    '**' : 2,
}

OPS_state = {
    'exp': 1,
    'pow2' : 1,
    'pow3' : 1,
    'sqrt' : 1,
    '-' : 1,
    '+' : 2,
    '*' : 2,
    '/' : 2,
    '**' : 2,
}

# Important metadata
metadata = {
            "Local_path": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results",
            "Exp_name": r"PropAmox50points",
            "Initial_p": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
            "Exp_T": np.linspace(607, 773, 10),
            "W_values": np.linspace(0, 8, 50),
            "chemicals": ["C3H6", "N2", "ACN", "NBP", "ACO", "OBP", "NH3", "O2", "H2O"]
            }

metadata_noisy = {
            "Local_path": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results",
            "Exp_name": r"PropAmox50pointsNoisy",
            "Initial_p": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
            "Exp_T": np.linspace(607, 773, 10),
            "W_values": np.linspace(0, 8, 50),
            "chemicals": ["C3H6", "N2", "ACN", "NBP", "ACO", "OBP", "NH3", "O2", "H2O"]
            }

metadata_noisy_less = {
            "Local_path": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results",
            "Exp_name": r"PropAmox25pointsNoisy",
            "Initial_p": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
            "Exp_T": np.linspace(607, 773, 10),
            "W_values": np.linspace(0, 8, 25),
            "chemicals": ["C3H6", "N2", "ACN", "NBP", "ACO", "OBP", "NH3", "O2", "H2O"]
            }

BMS_kwargs_mono =   {
                    "prior_folder": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS\Prior",
                    "chosen_output": 1,
                    "scaling": None,
                    "npar": 4,
                    "ops": OPS_mono
                    }

BMS_kwargs_state =   {
                    "prior_folder": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS\Prior",
                    "noutputs": 2,
                    "scaling": None,
                    "npar": 12,
                    "ops": OPS_state
                    }

test_conditions_dic = {"T_experiments": np.linspace(613, 763, 10),
                       "W_values": np.linspace(0,8,50),
                       "chemicals": ["ACN", "N2", "ACO"],
                       "integrator": "LSODA",
                       "p0": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
                        }

#! Try one. There are some weird things in the integration
# perform_exp = PropeneAmmoxidation(os.path.join(metadata["Local_path"], metadata["Exp_name"]))
# perform_exp.data_generation(initial_p=metadata["Initial_p"], T_experiments=metadata["Exp_T"],
#                             W_values=metadata["W_values"], integrator ="LSODA")
# perform_exp.profile_BMS_generation(chemicals =metadata["chemicals"], BMS_kwargs = BMS_kwargs_mono)
# perform_exp.data_state_BMS_generation(chemicals=metadata["chemicals"])
# perform_exp.run_state_BMS(data_folder = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results\PropAmox50points\BMS_state_space\nonnoisy",
#                           chemicals = metadata["chemicals"], BMS_kwargs=BMS_kwargs_state)
# perform_exp.calculate_test_experiment(BMS_models = "BMS_results3.xlsx", **test_conditions_dic)


#! Same with noise
# perform_exp2 = PropeneAmmoxidation(os.path.join(metadata_noisy["Local_path"], metadata_noisy["Exp_name"]), noise=True)
# perform_exp.data_generation(initial_p=metadata_noisy["Initial_p"], T_experiments=metadata_noisy["Exp_T"],
#                             W_values=metadata_noisy["W_values"], integrator ="LSODA")
# perform_exp.profile_BMS_generation(chemicals =metadata_noisy["chemicals"], BMS_kwargs = BMS_kwargs_mono)
# perform_exp.data_state_BMS_generation(chemicals=metadata_noisy["chemicals"])
# perform_exp.run_state_BMS(chemicals = metadata_noisy["chemicals"], BMS_kwargs=BMS_kwargs_state)

# Trying reading already in data
# perform_exp.data_state_BMS_generation(data = r"Local\Results\PropAmox50points\data\nonnoisy\data.xlsx",
#                                       bms_data_path=r"Local\Results\PropAmox50points\BMS_conc_profiles\nonnoisy\BMS_results.xlsx",
#                                       chemicals=["ACN", "ACO"])
# perform_exp2.calculate_test_experiment(BMS_models = "BMS_results.xlsx", **test_conditions_dic)

#! Now with less points but still noise
# perform_exp = PropeneAmmoxidation(os.path.join(metadata_noisy_less["Local_path"], metadata_noisy_less["Exp_name"]), noise=True)
# perform_exp.data_generation(initial_p=metadata_noisy_less["Initial_p"], T_experiments=metadata_noisy_less["Exp_T"],
#                             W_values=metadata_noisy_less["W_values"], integrator ="LSODA")
# perform_exp.profile_BMS_generation(chemicals =metadata_noisy_less["chemicals"], BMS_kwargs = BMS_kwargs_mono)
# perform_exp.data_state_BMS_generation(chemicals=metadata_noisy_less["chemicals"])
# perform_exp.run_state_BMS(chemicals = metadata_noisy_less["chemicals"], BMS_kwargs=BMS_kwargs_state)
# perform_exp.calculate_test_experiment(BMS_models = "BMS_results.xlsx", **test_conditions_dic)

#! Test integrators. LSODA works by changing the tolerance in the end
# try2 = PropeneAmmoxidation(os.path.join(metadata["Local_path"], "TestIntegration"))
# try2.data_generation(initial_p=metadata["Initial_p"], T_experiments=metadata["Exp_T"],
#                      W_values=np.linspace(0,10,500), integrator = "LSODA")

#! Various test
perform_exp = PropeneAmmoxidation(os.path.join(metadata["Local_path"], metadata["Exp_name"]), noise=False)
perform_exp.calculate_test_experiment(BMS_models = "BMS_results2.xlsx", **test_conditions_dic)