import numpy as np

#! Operations for the
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

#! Propene ammoxidation metadata
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

#! Methanol metadata
MeOH_metadata_nonnoisy_50 = {
            "Local_path": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results",
            "Exp_name": r"MeOH50points",
            "Initial_x": np.array([4, 0, 0, 82, 3, 11]), # %
            "Exp_T": np.linspace(450, 550, 10), 
            #"W_values": np.linspace(0, 34.8e-3, 50),
            "W_values": np.linspace(0, 0.25, 50),
            "chemicals": ["CO", "H2O", "CH3OH", "H2", "CO2", "Ar"]
            }

MeOHnn50_BMS_kwargs_mono =   {
                    "prior_folder": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS\Prior",
                    "chosen_output": 1,
                    "scaling": None,
                    "npar": 4,
                    "ops": OPS_mono
                    }

MeOHnn50_BMS_kwargs_state =   {
                    "prior_folder": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS\Prior",
                    "noutputs": 2,
                    "scaling": None,
                    "npar": 12,
                    "ops": OPS_state
                    }
