from scipy.integrate import odeint 
from data_generation import haber_bosch, ammonia_oxidation, propene_ammoxidation
import matplotlib.pyplot as plt
import numpy as np
from pyDOE import lhs
import pandas as pd
from sklearn.model_selection import train_test_split

class Experiment():
    def __init__(self, kinetic_function):
        self.kinetic_function = kinetic_function
        
    def sample_points(self, npoints:int, list_T:list, x_span:float, noise:float = False) -> tuple:
        """
        Samples a series of points using LHS for a number of isothermal temeperatures, for a given span of the 
        rate unit.
        
        Inputs:
            - npoints: Integer indicating the number of sampling points per experiment.
            - list_T: List of floats indicating the temperatures for which the experiments will take place.
            - x_span: Float that indicates the maximum advance of the reaction.
            - noise: If it exists, float that indicates the standard deviation of the normal noise introduced.
            
        Outputs:
            - true_results: Dictionary which has the elements of list_T as keys, and indicates the result of the 
            experiments for the whole integration space. Each element of the dictionary is a tuple, where in index 0
            the advance of reaction measurement, i.e., volume, amount of catalyst, etc., is indicated, and in index 1
            a numpy array of the concentration/pressure profile is shown.
            
            - sample_results: Dictionary which has the elements of list_T as keys and indicates the result of the experiments
            in the sampling points. Each element of the dictionary is a tuple, where in index 0
            the advance of reaction measurement, i.e., volume, amount of catalyst, etc., is indicated, and in index 1
            a numpy array of the concentration/pressure profile is shown.
            
            - noisy_results: Equivalent to sample_results, but the profile of integration results has noise added.
        """
        true_results = {i: self.perform_experiment(i) for i in list_T}
        np.random.seed(1)
        sample_x_points = lhs(1, npoints).flatten()
        sample_x_points = np.sort((sample_x_points*x_span))
        sample_indices = [0]
        for point in sample_x_points:
            closest_point = ((true_results[list(true_results.keys())[0]][0] - point)**2).argmin()
            if closest_point not in sample_indices:
                sample_indices.append(closest_point)
        sample_results = {} 
        sample_results_noisy = {}
        for experiment in true_results:
            W_sample = true_results[experiment][0][sample_indices]
            P_sample = true_results[experiment][1][sample_indices,:]
            if noise:
                P_sample_noisy = P_sample + np.random.normal(0, noise, size = P_sample.shape)
                sample_results_noisy[experiment] = (W_sample, P_sample_noisy)
            else:
                sample_results_noisy[experiment] = (W_sample, P_sample)
            sample_results[experiment] = (W_sample, P_sample)
        return true_results, sample_results, sample_results_noisy
         
        
class HaberBoschExperiment(Experiment):
    F0      = 29781.2          # kgmol/h
    T0      = 306.3 + 273.15   # K
    P0      = 203.957          # bar
    xH20    = 0.624            # bar
    xN20    = 0.183            # bar
    xNH30   = 0.136            # bar
    Lreactor = 2.13            # m
    Dreactor = 2               # m
    xInerts = 1 - xH20 - xN20 - xNH30
    p_cat  = 2200              # kgCat/m3
    R = 8.31                   # J/mol·K
    void_fraction = 0.33
    comp_list = ["Nitrogen", "Hydrogen", "Ammonia"]
    
    @staticmethod
    def ode_function(F, W, T, P, Finert, kinetic_function):
        FN2, FH2, FNH3 = F
        FT = FN2 + FH2 + FNH3 + Finert
        pN2  = P*FN2/FT 
        pH2  = P*FH2/FT 
        pNH3 = P*FNH3/FT 
        dN2dt  = - kinetic_function([pN2, pH2, pNH3], T)/2 
        dH2dt  = -3*kinetic_function([pN2, pH2, pNH3], T)/2
        dNH3dt = kinetic_function([pN2, pH2, pNH3], T)
        
        return [dN2dt, dH2dt, dNH3dt]

    def perform_experiment(self, T):
        Vreactor = HaberBoschExperiment.Lreactor * HaberBoschExperiment.Dreactor**2/4*np.pi
        Wreactor = Vreactor*HaberBoschExperiment.p_cat*(1-HaberBoschExperiment.void_fraction) # kg cat
        Finert = HaberBoschExperiment.P0*HaberBoschExperiment.xInerts
        W = np.linspace(0, Wreactor, 1000)
        F0 = [HaberBoschExperiment.F0*j for j in [HaberBoschExperiment.xN20,
                                                  HaberBoschExperiment.xH20,
                                                  HaberBoschExperiment.xNH30]]
        args = (T, HaberBoschExperiment.P0, Finert, self.kinetic_function)
        solve_ode = odeint(HaberBoschExperiment.ode_function, F0, W, args)
        return W, solve_ode

    def plot_experiment(self, results):
        plt.figure(dpi = 250)
        x_axis, y_axis = results
        for column in range(y_axis.shape[1]):
            plt.plot(x_axis, y_axis[:,column], label = HaberBoschExperiment.comp_list[column])
        plt.ylabel("F [kgmole/h]")
        plt.xlabel("W [kg Cat]")
        plt.legend()
        plt.show()

    def plot_multiple_experiment(self, results, rate_results, samples):
        plt.figure(dpi = 150)
        for n,figure in enumerate(HaberBoschExperiment.comp_list):
            plt.subplot(2,2, n+1)
            plt.title(f"Comp = {figure}")
            ax = plt.gca()
            for experiment in results:
                color = next(ax._get_lines.prop_cycler)['color']
                plt.plot(results[experiment][0], results[experiment][1][:,HaberBoschExperiment.comp_list.index(figure)], 
                         label = experiment, color = color)
                if samples:
                    plt.plot(samples[experiment][0], samples[experiment][1][:, HaberBoschExperiment.comp_list.index(figure)],
                             "o", color = color)
                plt.ylabel("F [kgmole/h]")
                plt.xlabel("W [kg Cat]")
            plt.legend()

        plt.subplot(2,2,4)
        ax = plt.gca()
        plt.title(f"Predicted vs True rate")
        true_rate, noisy_rate = rate_results 
        for experiment in true_rate:
            color = next(ax._get_lines.prop_cycler)['color']
            plt.plot(samples[experiment][0], true_rate[experiment], "--", color = color)
            plt.plot(samples[experiment][0][:-1], noisy_rate[experiment], "o", color = color)
        plt.xlabel("W [kg Cat]")
        plt.ylabel("$r_{NH_3} [\\frac{kgmole NH_3}{kg cat · h}]$")
        plt.show()
        
    def calculate_rate_samples(self, samples: dict, samples_noisy:dict , method = "forward"):
        """
        Converts the concentration units to partial pressures, and calculates the rate both with the true
        expression and a given method. For the forward method, calculates the rate at point "w" using the 
        partial pressures at point "w" and "w+1", calculating therefore the slope. Due to this, the dimension
        of true_rate and noisy_rate may differ.
        
        Some others methods that may be implemented are "backward", "forward-backward", and "polynomial".  
        """
        Finert = HaberBoschExperiment.P0*HaberBoschExperiment.xInerts
        true_rate  = {}
        noisy_rate = {}
        for experiment in samples:
            FT = samples[experiment][1].sum(axis = 1) + Finert
            x  = samples[experiment][1]/FT[:,None] 
            p  = x*HaberBoschExperiment.P0 
            true_rate[experiment] = self.kinetic_function(C = [p[:,0], p[:,1], p[:,2]], T = experiment)
        for experiment in samples_noisy:
            if method == "forward":
                noisy_rate[experiment] = np.diff(samples_noisy[experiment][1][:,2])/np.diff(samples_noisy[experiment][0])
        return true_rate, noisy_rate
    
    def generate_training_test(self, samples, rate, method = "forward", file_path = None, random_state = 1, train_size = 0.8):
        """
        Generates training and test data from samples and rates. For this case,
        the samples have to be in partial pressure, therefore they are changed. Then, the rates calculated
        from the calculate_rate_samples method are put together with the partial pressures and temperatures in a dataframe,
        which is written onto an .xlsx file if a file_path is given. 
        
        It repeats some functionality with calculate_rate_samples. Would be better to compress this functionality in only
        one method, to avoid having errors. 
        """
        dataframe = pd.DataFrame(columns = ["pN2 [bar]", "pH2 [bar]", "pNH3 [bar]", "T [K]", "rate [kgmole NH3/kgcat·h]"])
        Finert = HaberBoschExperiment.P0*HaberBoschExperiment.xInerts
        partial_pressures_list = []
        
        for experiment in samples:
            FT = samples[experiment][1].sum(axis = 1) + Finert 
            x = samples[experiment][1]/FT[:,None]
            p = x*HaberBoschExperiment.P0
            if method == "forward":
                input_data = p[:-1, :]
                T_vector   = np.repeat(experiment, input_data.shape[0]).reshape(-1,1)
                input_data = np.hstack((input_data, T_vector, rate[experiment][:input_data.shape[0]].reshape(-1,1)))
            partial_pressures_list.append(input_data)
        inputs = np.vstack(partial_pressures_list)
        dataframe[["pN2 [bar]", "pH2 [bar]", "pNH3 [bar]", "T [K]", "rate [kgmole NH3/kgcat·h]"]] = inputs
        train, test = train_test_split(dataframe, train_size=train_size, random_state=random_state)
        train = train.reset_index(drop = True)
        test  = test.reset_index(drop = True)
        if file_path:
            with pd.ExcelWriter(file_path) as writer:  
                train.to_excel(writer, sheet_name='train')
                test.to_excel(writer, sheet_name='test')
        return train,test
                                   


def main():
    HBE = HaberBoschExperiment(haber_bosch)
    true_results, sample_results, noisy_results = HBE.sample_points(npoints = 25, 
                                                    list_T = [450 + 10*i for i in range(25)], 
                                                    x_span = 4000,
                                                    noise = 20)
    true_rate, noisy_rate = HBE.calculate_rate_samples(sample_results, noisy_results)
    HBE.plot_multiple_experiment(true_results, rate_results = (true_rate, noisy_rate), samples = noisy_results)
    noisy_df = HBE.generate_training_test(noisy_results, noisy_rate, 
                                          file_path = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Data\Haber_Bosch\noisy_data.xlsx")
    true_df  = HBE.generate_training_test(sample_results, true_rate,
                                          file_path = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Data\Haber_Bosch\nonnoisy_data.xlsx")
if __name__ == "__main__":
    main()
