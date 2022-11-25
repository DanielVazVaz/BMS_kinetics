import numpy as np
import pandas as pd
from scipy.integrate import odeint
from pyDOE import lhs as py
import matplotlib.pyplot as plt

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
    
    
def propene_ammoxidation_kinetic(F, P:float,  T: float, Finert:float, BMS_ACN_eq = None):
    """Reference: https://www.sciencedirect.com/science/article/pii/S0021951716300239

    Args:
        F (np.array): np.array of the molar flows [mol/s]
        P (float): Total pressure [kPa]
        T (float): Temperature [K]
        Finert (float): Flow of inert [mol/s]

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
    if not BMS_ACN_eq:
        rACN    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*(1/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5)) + kACON*pACO*pNH3
    else:
        rACN    = BMS_ACN_eq(pC3H6, pN2, pACN, pNBP, pACO, pOBP, pNH3, pO2, pH2O, T)
    rNBP    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*((delta*pNH3 + epsilon*pO2**0.5)/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))
    rACO    = 0.94*kRLS*pC3H6*((1/(1 + alpha*pNH3))*(1/(1 + beta*pO2**0.5)) + (alpha*pNH3/(1 + alpha*pNH3))*(gamma*pH2O/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))) - kACON*pACO*pNH3
    rOBP    = 0.94*kRLS*pC3H6*(1/(1 + alpha*pNH3))*((beta*pO2**0.5)/(1 + beta*pO2**0.5)) + 0.06*kRLS*pC3H6
    rNH3    = -rACN - 7/3*rNBP - 2*rN2
    rO2     = -1.5*rACN - rACO - 7/3*rNBP - 10/3*rOBP - 3/2*rN2
    rH2O    = 3*rACN + rACO + 14/3*rNBP + 7/3*rOBP + 3*rN2
    
    return ([rC3H6, rN2, rACN, rNBP, rACO, rOBP, rNH3, rO2, rH2O],
            [pC3H6, pN2, pACN, pNBP, pACO, pOBP, pNH3, pO2, pH2O])

def ode_function(F, W: float, P: float, T:float, Finert: float):
    F[F < 0] = 0   # Check to avoid negative things
    return propene_ammoxidation_kinetic(F, P, T, Finert)[0]
    
def generate_experiment(F0, P:float, T: float):
    """Generation of an experiment

    Args:
        F0 (np.array): Initial vector of flows [mol/s]
        P (float): Isobaric pressure [kPa]
        T (float): Isothermal temperature [K]
    """    
    area_catalyst = 10          # m^2
    W   = np.linspace(0, area_catalyst, 500)
    Finert = calc_F_inert(F0)
    integration_results = odeint(ode_function, F0, W, args = (P, T, Finert))
    partial_pressures = integration_results/(integration_results.sum(axis = 1).reshape(-1,1) + Finert)*P
    return integration_results, W, partial_pressures, Finert

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

def plot_experiment(df, column_labels: list):
    """_summary_

    Args:
        df (_type_): _description_
        column_labels (list): _description_
    """    
    
    
    


    

