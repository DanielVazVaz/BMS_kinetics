from scipy.integrate import odeint 
from data_generation import haber_bosch, ammonia_oxidation, propene_ammoxidation
import matplotlib.pyplot as plt
import numpy as np
from pyDOE import lhs

class Experiment():
    def __init__(self, kinetic_function):
        self.kinetic_function = kinetic_function
         
        
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
        plt.show()
        
    def sample_points(self, npoints, list_T, x_span:int, noise = False):
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
    
    def calculate_rate_samples(self, samples, samples_noisy, method = "forward"):
        """
        Converts the concentration units to partial pressures, and the rate remains correct.
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
                
            
        
        
        
HBE = HaberBoschExperiment(haber_bosch)
true_results, sample_results, noisy_results = HBE.sample_points(npoints = 20, 
                                                list_T = [500, 550, 600, 650, 700], 
                                                x_span = 9000,
                                                noise = 50)
true_rate, noisy_rate = HBE.calculate_rate_samples(sample_results, noisy_results)
HBE.plot_multiple_experiment(true_results, rate_results = (true_rate, noisy_rate), samples = noisy_results)

