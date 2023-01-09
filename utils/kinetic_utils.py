import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, odeint
import json
import sympy as sym
import matplotlib.pyplot as plt
from .BMS_utils import BMS_instance
import os
from itertools import product

class PropeneAmmoxidation:
    """Dummy class to maintain these functions together
    """
    def __init__(self, dump_folder:str, noise = False):
        self.comp_order = ["C3H6", "N2", "ACN", "NBP", "ACO", "OBP", "NH3", "O2", "H2O"]
        self.dump_folder = dump_folder
        self.noise = noise 
        
    def data_generation(self, initial_p, T_experiments, W_values = None, file_name = r"data.xlsx", integrator = "LSODA"):
        """Generates the data
        """
        if file_name[-5:] != ".xlsx":
            file_name += ".xlsx"
        # Create the folder to hold the data in case it does not exist
        data_folder = os.path.join(self.dump_folder, "data")
        if self.noise:
            folder_to_save = os.path.join(data_folder, "noisy")
        else:
            folder_to_save = os.path.join(data_folder, "nonnoisy")
        os.makedirs(folder_to_save, exist_ok=True)
        
        # Calculate the initial flow rate
        P0 = 120    # kPa
        T0 = 298    # K
        v0 = 100    # mL/min
        R  = 0.082  # atm·L/mol·K
        p0 = initial_p
        F0 = p0/101.3*v0/1000/(R*T0)/60
        
        labels = ["$C_3H_6$", "$N_2$", "$ACN$", "$NBP$", "$ACO$", "$OBP$", "$NH_3$", "$O_2$", "$H_2O$"]
        df_dic = {}
        for n,T in enumerate(T_experiments):
            profile, W, partial, Finert = PropeneAmmoxidation.generate_experiment(F0, P = P0, T = T, W_values=W_values, integrator=integrator)
            df_dic[T] = PropeneAmmoxidation.generate_df(labels, T, partial, W, profile, Finert)
            multi_index = pd.MultiIndex.from_tuples(product([f"Exp{n+1}"], df_dic[T].index))
            df_dic[T].set_index(multi_index, inplace=True)
        df_complete = pd.concat([df_dic[i] for i in df_dic])
        if self.noise:
            chemicals = self.comp_order
            chemical_flows = ["F" + i + " [mol/s]" for i in chemicals]
            P_flows = ["p" + i + " [kPa]" for i in chemicals]
            np.random.seed(1)
            noise_value = np.random.normal(loc = 1, scale = 0.05, size = df_complete.loc[:, chemical_flows].shape)
            np.random.seed(None)
            df_complete.loc[:, chemical_flows] = df_complete.loc[:, chemical_flows]*noise_value
            df_complete.loc[:, P_flows] = df_complete.loc[:, P_flows]*noise_value
        self.data_exp = df_complete
        df_complete.to_excel(os.path.join(folder_to_save, file_name))
        
    def profile_BMS_generation(self, data = None, chemicals = None, max_conv = 100, BMS_kwargs = None):
        """Generates the monodimensional BMS that will afterwards be differentiated.

        Args:
            data (str, optional): Data of the experiments. Defaults to None.
            chemicals (list, optional): List of chemicals to do the BMS. Defaults to None.
            max_conv (float, optional): Maximum conversion allowed of C3H6. Defaults to 100.
            BMS_kwargs (dict, optional): Arguments to the BMS instance. Defaults to None.
        """        
        if not data:
            data = self.data_exp
        else:
            data = pd.read_excel(data, index_col = [0,1]) 
        
        # Create the folder to hold the BMS excel equations in case it does not exist
        data_folder = os.path.join(self.dump_folder, "BMS_conc_profiles")
        if self.noise:
            folder_to_save = os.path.join(data_folder, "noisy")
        else:
            folder_to_save = os.path.join(data_folder, "nonnoisy")
        os.makedirs(folder_to_save, exist_ok=True)
        
        # Filter regarding the conversion
        self.data_filtered = data[data["convC3H6 [%]"] <= max_conv]
        
        # Check experiments
        exp_data_folder = folder_to_save
        os.makedirs(exp_data_folder, exist_ok=True)
        BMS_instances = {}
        for experiment in self.data_filtered.index.get_level_values(0).unique():
            print(f"[+] Experiment: {experiment}")
            exp_data = self.data_filtered.loc[experiment]
            for chemical in chemicals:
                chemical_flow = "F" + chemical + " [mol/s]"
                exp_df = exp_data[["W [m2]", chemical_flow]]
                exp_file_name = f"{experiment.replace(' ','')}_{chemical}.xlsx"
                data_BMS_path = os.path.join(exp_data_folder, exp_file_name)
                with pd.ExcelWriter(data_BMS_path) as writer:
                    exp_df.to_excel(writer, sheet_name=  "train")
                    exp_df.to_excel(writer, sheet_name = "test")
                BMS_instances[experiment, chemical] = BMS_instance(data_path = data_BMS_path, **BMS_kwargs)
                BMS_instances[experiment, chemical].run_BMS(mcmcsteps = 500)
                ntries = 0
                R2_train = BMS_instances[experiment,chemical].calculate_rsquare()["R2_train"]
                while (ntries <= 3) and (R2_train <= 0.925):
                    print(f"Recalculating. Actual R2 = {R2_train}")
                    BMS_instances[experiment, chemical].run_BMS(mcmcsteps = 500)
                    R2_train = BMS_instances[experiment,chemical].calculate_rsquare()["R2_train"]
                    ntries += 1
                    
        self.BMS_instances = BMS_instances
        
        # We create the dataframe structure to save the BMS equations
        dic_BMS_eq = {}
        x1 = sym.symbols("x1")
        for chemical in chemicals:
            dic_BMS_eq[chemical] = pd.DataFrame(index = self.data_filtered.index.get_level_values(0).unique(),
                                     columns = ["Model", "ParVal", "R2"])
            for experiment in dic_BMS_eq[chemical].index:
                dic_BMS_eq[chemical].loc[experiment, "Model"] = str(BMS_instances[experiment, chemical].mdl_model)
                dic_BMS_eq[chemical].loc[experiment, "ParVal"] = str(BMS_instances[experiment, chemical].mdl_model.par_values["d0"])
                dic_BMS_eq[chemical].loc[experiment, "R2"] = BMS_instances[experiment, chemical].calculate_rsquare()["R2_train"]
                # Differential calculation:
                dequation = sym.diff(dic_BMS_eq[chemical].loc[experiment, "Model"], x1)
                dequation = dequation.subs(BMS_instances[experiment, chemical].mdl_model.par_values["d0"])
                dic_BMS_eq[chemical].loc[experiment, "rate"] = str(dequation)
                
        with pd.ExcelWriter(os.path.join(exp_data_folder,"BMS_results.xlsx")) as writer:
            for chemical in dic_BMS_eq.keys():
                dic_BMS_eq[chemical].to_excel(writer, sheet_name = chemical)
        self.dic_BMS_eq = dic_BMS_eq
                
    def data_state_BMS_generation(self, data = None, bms_data_path = None,
                                  chemicals = None, max_conv = 100):
        """Generates the data for the BMS in the state space

        Args:
            data (str, optional): Data. Defaults to None.
            max_conv (float, optional): Maximum conversion of C3H6. Defaults to 100.
        """        
        if not data:
            data = self.data_exp
        else:
            data = pd.read_excel(data, index_col = [0,1]) 
            
        if not bms_data_path:
            bms_data = self.dic_BMS_eq
        else:
            bms_data = {}
            for chemical in chemicals:
                bms_data[chemical] = pd.read_excel(bms_data_path, sheet_name=chemical, index_col=0)
            self.dic_BMS_eq = bms_data 

        # Create the folder to hold the BMS data in case it does not exist
        data_folder = os.path.join(self.dump_folder, "BMS_state_space")
        if self.noise:
            folder_to_save = os.path.join(data_folder, "noisy")
        else:
            folder_to_save = os.path.join(data_folder, "nonnoisy")
        os.makedirs(folder_to_save, exist_ok=True)
        
        # Filter regarding the conversion
        self.data_filtered_2 = data[data["convC3H6 [%]"] <= max_conv]
        self.data_valid = self.data_filtered_2.copy()
        
        # Data for reading the bms equations
        x1 = sym.symbols("x1")
        
        # Iterate through experiments and chemicals to generate true and calculated rate
        matching_output_chemical = {i:j for i,j in zip(self.comp_order, range(9))}
        P0 = 120
        
        for row in self.data_valid.index:
            F = self.data_valid.loc[row, "FC3H6 [mol/s]":"FH2O [mol/s]"].to_numpy()
            T = self.data_valid.loc[row, "T [K]"]
            Finert = self.data_valid.loc[row, "Finert [mol/s]"]
            for chemical in chemicals:
                self.data_valid.loc[row, f"True r{chemical} [mol/m2s]"] = PropeneAmmoxidation.propene_ammoxidation_kinetic(F = F, T = T, Finert = Finert, 
                                                                                                                           P = P0)[0][matching_output_chemical[chemical]]
                diff1 = bms_data[chemical].loc[row[0], "rate"]
                bms_diff_expr = sym.lambdify(args = x1, expr = sym.sympify(diff1))
                self.data_valid.loc[row, f"Calc r{chemical} [mol/m2s]"] = bms_diff_expr(self.data_valid.loc[row, "W [m2]"])
                
        # Remove inf, -inf, and nan rows.
        self.data_valid = self.data_valid.replace([np.inf, -np.inf], np.nan).dropna()
        
        # Generate the excel files necessary for the BMS stuff 
        self.data_training_dic = {}   
        for chemical in chemicals:
            data_training = self.data_valid.loc[:, "pC3H6 [kPa]":"T [K]"].copy()
            data_training["Extra1"] = 0.0 # We do not have prior for 10 variables
            data_training["Extra2"] = 0.0 # We do not have prior for 10 variables
            data_training["BMS"]    = self.data_valid.loc[:, f"Calc r{chemical} [mol/m2s]"]
            data_training["True"]   = self.data_valid.loc[:, f"True r{chemical} [mol/m2s]"]
            data_training.index = [(i,j) for (i,j) in data_training.index]
            self.data_training_dic[chemical] = data_training 
            with pd.ExcelWriter(os.path.join(folder_to_save, f"BMS_train_{chemical}.xlsx")) as writer:
                data_training.to_excel(writer, sheet_name = "train")
                data_training.to_excel(writer, sheet_name = "test")
    
    def run_state_BMS(self, data_folder = None, chemicals = None, BMS_kwargs = None, mcmcsteps = 5000):
        """Generates the BMS for the state variables for the designed chemicals

        Args:
            data_folder (_type_, optional): Path to the folder where the training excels are. Defaults to None.
            chemicals (_type_, optional): List of chemicals. Defaults to None.
            BMS_kwargs (_type_, optional): Keyword arguments to the BMS creation. Defaults to None.
            mcmcsteps (int, optional): Number of mcmcsteps. Defaults to 5000.
        """        

        # Folder with all the BMS excel data
        if not data_folder:
            data_folder = os.path.join(self.dump_folder, "BMS_state_space")
            if self.noise:
                data_folder = os.path.join(data_folder, "noisy")
            else:
                data_folder = os.path.join(data_folder, "nonnoisy")
                
        # Folder to create to save the result of the BMS run. If it ever runs completely.
        save_folder = os.path.join(self.dump_folder, "BMS_state_results")
        if self.noise:
            save_folder = os.path.join(save_folder, "noisy")
        else:
            save_folder = os.path.join(save_folder, "nonnoisy")
        os.makedirs(save_folder, exist_ok=True)
        
        
        
        # Loop to create the BMS things
        self.BMS_state_dic = {}
        self.df_dic = {}
        for chemical in chemicals:
            excel_file = os.path.join(data_folder, f"BMS_train_{chemical}.xlsx")
            self.BMS_state_dic[chemical, "True"] = BMS_instance(data_path = excel_file, chosen_output = 2, **BMS_kwargs)
            self.BMS_state_dic[chemical, "True"].run_BMS(mcmcsteps=mcmcsteps)
            self.BMS_state_dic[chemical, "Calc"] = BMS_instance(data_path = excel_file, chosen_output = 1, **BMS_kwargs)
            self.BMS_state_dic[chemical, "Calc"].run_BMS(mcmcsteps=mcmcsteps)
            
            self.df_dic[chemical] = pd.DataFrame(index = ["True", "Calc"], columns = ["Model", "ParValues", "R2"])
            for model in self.df_dic[chemical].index:
                self.df_dic[chemical].loc[model, "Model"] = str(self.BMS_state_dic[chemical,model].mdl_model)
                self.df_dic[chemical].loc[model, "ParValues"] = str(self.BMS_state_dic[chemical,model].mdl_model.par_values["d0"])
                self.df_dic[chemical].loc[model, "R2"] = self.BMS_state_dic[chemical,model].calculate_rsquare()["R2_train"]
                self.df_dic[chemical].loc[model, "Time [min]"] = self.BMS_state_dic[chemical,model].training_time
        with pd.ExcelWriter(os.path.join(save_folder, f"BMS_results.xlsx")) as writer:
            for chemical in chemicals:
                self.df_dic[chemical].to_excel(writer, sheet_name = chemical)

    def calculate_test_experiment(self, BMS_models = "BMS_results.xlsx", chemicals = None, integrator = "LSODA",
                             T_experiments = None, W_values = None, p0 = None):
        """Plots a test experiment

        Args:
            BMS_models (path, optional): Path to excel file with the BMS models saved. Defaults to None.
            chemicals (_type_, optional): List of chemicals to represent. Defaults to None.
            integrator (str, optional): Integrator for the solve_ivp function. Defaults to "LSODA".
            T_experiments (list): Experiments to be done
            W_values: Values where the integrator will answer something
        """        
        
        # We set up the path of the xlsx
        BMS_models_folder = os.path.join(self.dump_folder, "BMS_state_results", "noisy" if self.noise else "nonnoisy")
        BMS_models = os.path.join(BMS_models_folder, BMS_models)
        
        matching_output_chemical = {i:j for i,j in zip(self.comp_order, range(9))}
        # We read the True and calc function of the excel file and numpify it
        self.BMS_True_eq = {}
        self.BMS_Calc_eq = {}
        for chemical in chemicals:
            df_to_read = pd.read_excel(BMS_models, sheet_name=chemical, index_col=0)
            self.BMS_True_eq[chemical] = sym.sympify(df_to_read.loc["True", "Model"]).subs(json.loads(df_to_read.loc["True", "ParValues"].replace("'", '"')))
            self.BMS_True_eq[chemical] = self.BMS_True_eq[chemical].subs({"x11": 0, "x12": 0})
            self.BMS_Calc_eq[chemical] = sym.sympify(df_to_read.loc["Calc", "Model"]).subs(json.loads(df_to_read.loc["Calc", "ParValues"].replace("'", '"')))
            self.BMS_Calc_eq[chemical] = self.BMS_Calc_eq[chemical].subs({"x11": 0, "x12": 0})
            self.BMS_True_eq[chemical] = sym.lambdify(["x" + str(i+1) for i in range(10)], self.BMS_True_eq[chemical])
            self.BMS_Calc_eq[chemical] = sym.lambdify(["x" + str(i+1) for i in range(10)], self.BMS_Calc_eq[chemical])

        # Initial conditions
        P0 = 120    # kPa
        T0 = 298    # K
        v0 = 100    # mL/min
        R  = 0.082  # atm·L/mol·K
        labels = ["$C_3H_6$", "$N_2$", "$ACN$", "$NBP$", "$ACO$", "$OBP$", "$NH_3$", "$O_2$", "$H_2O$"]
        #p0 = np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5])
        p0 = p0
        F0 = p0/101.3*v0/1000/(R*T0)/60
        
        self.test_experiments = {}
        self.match_true_eq = [(matching_output_chemical[i],self.BMS_True_eq[i]) for i in chemicals]
        self.match_calc_eq = [(matching_output_chemical[i],self.BMS_Calc_eq[i]) for i in chemicals]
        self.match_bms_eq = {"Real": None, "True": self.match_true_eq, "Calc": self.match_calc_eq}
        for value in ["Real", "True", "Calc"]:
            df_dic = {}
            for n,T in enumerate(T_experiments):
                # RLoop for the values
                profile, W, partial, Finert = PropeneAmmoxidation.generate_experiment(F0, P = P0, T = T, W_values=W_values, integrator=integrator, BMS_eq=self.match_bms_eq[value])
                df_dic[T] = PropeneAmmoxidation.generate_df(labels, T, partial, W, profile, Finert)
                multi_index = pd.MultiIndex.from_tuples(product([f"Exp{n+1}"], df_dic[T].index))
                df_dic[T].set_index(multi_index, inplace = True)
                if df_dic[T].isna().sum().sum() > 0:
                    print(f"Problem in integration in BMS {value} at experiment {df_dic[T].index[0][0]}")
                    
            self.test_experiments[value] = pd.concat([df_dic[i] for i in df_dic])
                
    
    def plot_test_experiment(self, experiment: str, chemicals:list):
        chemical_flows = ["F" + i + " [mol/s]" for i in chemicals]
        plt.figure(dpi = 300)
        for chemical, chemical_flow in zip(chemicals, chemical_flows):
            W_true = self.test_experiments["Real"].loc[experiment, "W [m2]"]
            F_true = self.test_experiments["Real"].loc[experiment, chemical_flow]
            plt.plot(W_true, F_true, label = chemical, ls = "-")
            for value in ["True", "Calc"]:
                F = self.test_experiments[value].loc[experiment, chemical_flow]
                notna_values = F.notna()
                F = F[notna_values]
                W = self.test_experiments[value].loc[experiment, "W [m2]"]
                W = W[notna_values]
                SST = sum((F - F_true[notna_values])**2)
                SSR = sum((F_true[notna_values] - F_true[notna_values].mean())**2)
                try:
                    R2  = 1 - SST/SSR
                except:
                    R2 = "nan"
                plt.plot(W, F, label = f"{value} ({R2:.4f})", ls = "--")
        plt.xlabel("W [$m^2$]")
        plt.ylabel("F [mol/s]")
        plt.legend()
        plt.title(experiment)
        plt.show()
    @staticmethod
    def calc_F_inert(F):
        """Calculates the amount of inert in the reaction. The assumption is that the volumetric
        flow at 273 K is equal to 100 mL/min.

        Args:
            F (np.array): np.array of the molar flows [mol/s]
        """    
        FT = F.sum()     # mol/s
        R = 0.082        # atm·L/mol·K
        v = 100          # mL/min
        T = 298          # K
        P = 120          # kPa
        Finert = P/101.3*v/1000/(R*T)/60 - FT    # mol/s
        assert Finert > 0
        return Finert
        
    @staticmethod
    def propene_ammoxidation_kinetic(F, P:float,  T: float, Finert:float, BMS_eq = None):
        """Reference: https://www.sciencedirect.com/science/article/pii/S0021951716300239

        Args:
            F (np.array): np.array of the molar flows [mol/s]
            P (float): Total pressure [kPa]
            T (float): Temperature [K]
            Finert (float): Flow of inert [mol/s]
            BMS_eq (listof tuples): Equation in case it comes from the BMS. e.g. BMS_eq = [(0, rC3G6eq), (1, rN2eq)]

        Returns:
            tuple: value of the different rates [mol/s·m^2]
        """
        F[F<0] = 0
        FC3H6, FN2, FACN, FNBP, FACO, FOBP, FNH3, FO2, FH2O = F                         # mol/s
        pC3H6, pN2, pACN, pNBP, pACO, pOBP, pNH3, pO2, pH2O = F/(F.sum() + Finert)*P    # kPa
        
        R       = 1.9858775e-3                  # kcal/mol·K
        kRLS    = 2.3*np.exp(-22.4/(R*T))       # mol/s/m2/kPa
        kACON   = 6e-4*np.exp(-7/(R*T))         # mol/s/m2/kPa2
        kN2     = 9e5*np.exp(-35/(R*T))         # mol/s/m2/kPa
        kH2O    = 5e1*np.exp(-2/(R*T))          # mol/s/m2/kPa
        alpha   = 3e-5*np.exp(15/(R*T))         # mol/s/m2/kPa
        beta    = 2e-3*np.exp(5/(R*T))          # mol/s/m2/kPa
        gamma   = 2*np.exp(-10/(R*T))           # mol/s/m2/kPa
        delta   = 1e-5*np.exp(10/(R*T))         # mol/s/m2/kPa
        epsilon = 2e-3*np.exp(5/(R*T))          # mol/s/m2/kPa
        
        # Rates are in mol/m2·s
        rC3H6   = -kRLS*pC3H6
        rN2     = kN2*pNH3/(1 + kH2O*pH2O)
        rACN    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*(1/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5)) + kACON*pACO*pNH3
        rNBP    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*((delta*pNH3 + epsilon*pO2**0.5)/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))
        rACO    = 0.94*kRLS*pC3H6*((1/(1 + alpha*pNH3))*(1/(1 + beta*pO2**0.5)) + (alpha*pNH3/(1 + alpha*pNH3))*(gamma*pH2O/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))) - kACON*pACO*pNH3
        rOBP    = 0.94*kRLS*pC3H6*(1/(1 + alpha*pNH3))*((beta*pO2**0.5)/(1 + beta*pO2**0.5)) + 0.06*kRLS*pC3H6
        rNH3    = -rACN - 7/3*rNBP - 2*rN2
        rO2     = -1.5*rACN - rACO - 7/3*rNBP - 10/3*rOBP - 3/2*rN2
        rH2O    = 3*rACN + rACO + 14/3*rNBP + 7/3*rOBP + 3*rN2
        
        rate_list = [rC3H6, rN2, rACN, rNBP, rACO, rOBP, rNH3, rO2, rH2O]
        
        if BMS_eq:
            for eq in BMS_eq:
                rate_list[eq[0]] = eq[1](pC3H6, pN2, pACN, pNBP, pACO, pOBP, pNH3, pO2, pH2O, T)
        
        return (rate_list,
                [pC3H6, pN2, pACN, pNBP, pACO, pOBP, pNH3, pO2, pH2O])
        
    @staticmethod 
    def ode_function(W, F, P: float, T:float, Finert: float, BMS_eq = None):
        F[F < 0] = 0   # Check to avoid negative things
        return PropeneAmmoxidation.propene_ammoxidation_kinetic(F, P, T, Finert, BMS_eq)[0]
        
    @staticmethod
    def generate_experiment(F0, P:float, T: float, W_values = None, integrator = "LSODA", BMS_eq = None):
        """Generation of an experiment

        Args:
            F0 (np.array): Initial vector of flows [mol/s]
            P (float): Isobaric pressure [kPa]
            T (float): Isothermal temperature [K]
        """    
        area_catalyst = 10          # m^2
        if W_values is None:
            W   = np.linspace(0, area_catalyst, 500)
        else:
            W = W_values
        Finert = PropeneAmmoxidation.calc_F_inert(F0)
        if integrator == "ODEINT":
            integration_results = odeint(PropeneAmmoxidation.ode_function, F0, W, args = (P,T,Finert, BMS_eq), tfirst=True)
            y_integration = integration_results
            partial_pressures = integration_results/(integration_results.sum(axis = 1).reshape(-1,1) + Finert)*P
        else:
            integration_results = solve_ivp(PropeneAmmoxidation.ode_function, 
                                            (0, area_catalyst), F0, t_eval = W, args = (P, T, Finert, BMS_eq),
                                            method=integrator, rtol=1.49012e-8, atol=1.49012e-8)
            y_integration = integration_results.y.T
            partial_pressures = y_integration/(y_integration.sum(axis = 1).reshape(-1,1) + Finert)*P
        return y_integration, W, partial_pressures, Finert

    @staticmethod
    def generate_df(latex_labels: str, T: float, partial_pressures: any, W:any, F: any, Finert: float) -> any:
        """Generate the dataframes for a experiment

        Args:
            latex_labels (str): Latex label for components, i.e., "$CO_2$"
            T (float): Temperature of the experiment [K]
            partial_pressures (array): Partial pressures of the components [kPa]
            W (array): Vector of catalyst surface [m2]
            F (array): Molar flows of the profiles [mol/s]
            Finert (float): Molar flow of the inert in this experiment [mol/s]

        Returns:
            any: Pandas dataframe of a given experiment
        """
        df_plabels = ["p" + i[1:-1].replace("_","") + " [kPa]" for i in latex_labels]
        df_flabels = ["F" + i[1:-1].replace("_","") + " [mol/s]" for i in latex_labels]
        df_labels = df_plabels + ["T [K]", "W [m2]"] + df_flabels
        df = pd.DataFrame(index = range(partial_pressures.shape[0]), columns = df_labels)
        df.loc[:, df_plabels] = partial_pressures
        df.loc[:, "T [K]"] = T
        df.loc[:, "W [m2]"] = W 
        df.loc[:, df_flabels] = F
        df.loc[:, "Finert [mol/s]"] = Finert
        # We check the conversion of C3H6
        df.loc[:, "convC3H6 [%]"] = (df.loc[0, "FC3H6 [mol/s]"] - df.loc[:, "FC3H6 [mol/s]"])/(df.loc[0, "FC3H6 [mol/s]"])*100
        return df 
    
    @staticmethod
    def plot_experiment(data, chemicals:list):
        chemical_flows = ["F" + i + " [mol/s]" for i in chemicals]
        W = data.loc[:, "W [m2]"]
        plt.figure()
        for chemical,chemical_flow in zip(chemicals,chemical_flows):
            plt.plot(W, data.loc[:, chemical_flow], label = chemical)
        plt.xlabel("W [$m^2$]")
        plt.ylabel("F [$mol/s$]")
        plt.legend(loc = "best")
        plt.show()
        

        
            

            
            
            
        
