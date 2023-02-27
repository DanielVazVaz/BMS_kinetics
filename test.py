import sys
BMS_FOLDER = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS"
sys.path.append(BMS_FOLDER)
sys.path.append(BMS_FOLDER + "\\" + "Prior")
from utils.kinetic_utils import PropeneAmmoxidation, MethanolProduction
from metadata_file import (OPS_mono, OPS_state, MeOH_metadata_nonnoisy_50, 
                           MeOHnn50_BMS_kwargs_mono,MeOHnn50_BMS_kwargs_state)
import numpy as np
import os

# Important metadata


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
#perform_exp = PropeneAmmoxidation(os.path.join(metadata["Local_path"], metadata["Exp_name"]), noise=False)
#perform_exp.calculate_test_experiment(BMS_models = "BMS_results2.xlsx", **test_conditions_dic)

#!Methanol stuff
MeOH50nn = MethanolProduction(os.path.join(MeOH_metadata_nonnoisy_50["Local_path"], MeOH_metadata_nonnoisy_50["Exp_name"]))
# MeOH50nn.data_generation(initial_x=MeOH_metadata_nonnoisy_50["Initial_x"], T_experiments=MeOH_metadata_nonnoisy_50["Exp_T"],
#                             W_values=MeOH_metadata_nonnoisy_50["W_values"], integrator ="LSODA")
# MeOH50nn.plot_experiment(MeOH50nn.data_exp.loc["Exp5"], ["CO", "CO2", "H2O", "CH3OH"])
# MeOH50nn.profile_BMS_generation(chemicals =["CO", "CH3OH"], BMS_kwargs = MeOHnn50_BMS_kwargs_mono)
# MeOH50nn.data_state_BMS_generation(data = r"Local\Results\MeOH50points\data\nonnoisy\data.xlsx",
#                                   bms_data_path=r"Local\Results\MeOH50points\BMS_conc_profiles\nonnoisy\BMS_results.xlsx",
#                                   chemicals = ["CO", "CH3OH"])
# MeOH50nn.run_state_BMS(data_folder = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results\MeOH50points\BMS_state_space\nonnoisy",
#                        chemicals = ["CO", "CH3OH"], BMS_kwargs=MeOHnn50_BMS_kwargs_state)
MeOH50nn.calculate_test_experiment("BMS_results.xlsx", chemicals = ["CO", "CH3OH"], T_experiments=MeOH_metadata_nonnoisy_50["Exp_T"][:1], 
                                   W_values = MeOH_metadata_nonnoisy_50["W_values"], x0 = MeOH_metadata_nonnoisy_50["Initial_x"])

# MeOH50n = MethanolProduction(os.path.join(MeOH_metadata_nonnoisy_50["Local_path"],
#                                           MeOH_metadata_nonnoisy_50["Exp_name"]),
#                              noise = True)
# MeOH50n.data_generation(initial_x=MeOH_metadata_nonnoisy_50["Initial_x"], T_experiments=MeOH_metadata_nonnoisy_50["Exp_T"],
#                             W_values=MeOH_metadata_nonnoisy_50["W_values"], integrator ="LSODA")
# MeOH50n.profile_BMS_generation(chemicals =["CO", "CH3OH"], BMS_kwargs = MeOHnn50_BMS_kwargs_mono)
# MeOH50n.data_state_BMS_generation(chemicals = ["CO", "CH3OH"])
# MeOH50n.run_state_BMS(chemicals = ["CO", "CH3OH"], BMS_kwargs=MeOHnn50_BMS_kwargs_state)
