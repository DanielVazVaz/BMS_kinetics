#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from pyDOE import lhs
import pandas as pd
from sklearn.model_selection import train_test_split

def generate_excel(data: pd.DataFrame, folder: str, random_state = 1, train_size = 0.8):
    """
    Generates training and test sheets from the lhs data design.
    """
    train,test = train_test_split(data, random_state = random_state, train_size = train_size)
    train = train.reset_index(drop = True)
    test  = test.reset_index(drop = True)
    with pd.ExcelWriter(folder + "\\" + r"data.xlsx") as writer:  
        train.to_excel(writer, sheet_name='train')
        test.to_excel(writer, sheet_name='test')
    return train,test

def haber_bosch(C: list, T: float):
    """
    Reaction rate for the Haber-Bosch process in kgmol NH3/kgcat·h
    N2 + 3H2 <-> 2NH3
    
    Reference: https://www.sciencedirect.com/science/article/pii/S0098135408000392#app1
    """
    f = 4.75
    R = 8.31            # J/mol·K
    pcat = 2200         # kg/m3
    pN2, pH2, pNH3 = C  # bar
    k1 = 1.79e4*np.exp(-87090/(R*T))
    k2 = 2.75e16*np.exp(-198464/(R*T))
    # rNH3 = 2*f/pcat*(k1*(pN2*pH2**1.5)/(pNH3) - k2*pNH3/pH2**1.5) # This is reported in the paper
    rNH3 = 2*f*(k1*(pN2*pH2**1.5)/(pNH3) - k2*pNH3/pH2**1.5)   # For some reason, this is the one they use. I'm missing something here.  
    return rNH3

def haber_bosch_experiment(nsamples = 625):
    """
    Develops the initial conditions, the variation, and creates training
    and test data
    
    Reference: https://www.sciencedirect.com/science/article/pii/S0098135408000392#app1
    """
    T0    = 306.3 + 273.15   # K
    P0    = 203.957          # bar
    pH20  = P0*0.624         # bar
    pN20  = P0*0.183         # bar
    pNH30 = P0*0.136         # bar
    lower = [P0*0.527,P0*0.149,0.136*P0,T0]  # H2, N2, NH3, T
    upper = [P0*0.624,P0*0.183,0.262*P0,452.2+273.15]
    np.random.seed(1)
    lhs_data = lhs(4, nsamples)
    for column in range(lhs_data.shape[1]):
        lhs_data[:,column] = lhs_data[:,column]*(upper[column] - lower[column]) + lower[column]
    data = pd.DataFrame(lhs_data)
    data.columns = ["pH2 [bar]", "pN2 [bar]", "pNH3 [bar]", "T [K]"]
    data["rate [kgmole NH3/kgcat·h"] = haber_bosch([data["pN2 [bar]"],
                                                    data["pH2 [bar]"],
                                                    data["pNH3 [bar]"]],
                                                   data["T [K]"])
    train,test = generate_excel(data, folder = r"Local\Data\Haber_Bosch")
    return train,test

def ammonia_oxidation(C, T):
    """
    Reaction rate for the oxidation of amonia to NO and N2 in mol cm^-2 s^-1
    
    2 NH3 + 3/2 O2 -> N2 + 3 H2O
    2 NH3 + 5/2 O2 -> 2 NO + 3 H2O
    
    Reference: https://www.sciencedirect.com/science/article/pii/0021951775902493?via%3Dihub
    """
    pNH3, pO2, pNO = C              # Torr
    R = 8.31                        # J/mol·K
    
    rN2 = (
            (4.2e-9*np.exp(26_000/(R*T))*pNO*pNH3)/(1 + 5.4e-2*np.exp(10_900/(R*T))*pO2 + 4.5e-6*np.exp(20_900/(R*T))*pNH3)**2
            +
            (2.2e-6*np.exp(-570/(R*T))*pNH3)/(1 + 5.4e-2*np.exp(10_900/(R*T))*pO2 + 4.5e-6*np.exp(20_900/(R*T))*pNH3)
          )
    
    rNO = 3.4e-8*np.exp(21_700/(R*T))*pNH3*(pO2**(1/2))/((1+8e-2*np.exp(4400/(R*T))*pO2**(1/2)) * (1 + 1.6e-3*np.exp(25_500/(R*T))*pNH3))
    return rN2, rNO

def ammonia_oxidation_experiment(nsamples = 625):
    """
    Develops the initial condition, variations, and things that must be taken into account.
    
    The total pressure must be between 0.1 and 1 Torr. We will use an LHS design  and maintain this with random.
    
    The temperature will vary from 1000 K to 1273 K.
    
    Reference: https://www.sciencedirect.com/science/article/pii/0021951775902493?via%3Dihub
    """
    lower = [0.02, 0.02, 0.02, 1173]   # pNH3, pO2, pNO, T
    upper = [1, 1, 1, 1273]         # pNH3, PO2, pNO, T
    np.random.seed(1)
    lhs_data = lhs(4, nsamples)
    for column in range(lhs_data.shape[1]):
        lhs_data[:,column] = lhs_data[:,column]*(upper[column] - lower[column]) + lower[column]
    for row in range(lhs_data.shape[0]):
        total_P = lhs_data[row, :3].sum()
        if total_P > 1:
            lhs_data[row, :3] /= total_P 
        elif total_P < 0.1:
            lhs_data[row, :3] /= total_P*0.1
    data = pd.DataFrame(lhs_data)
    data.columns = ["pNH3 [Torr]", "pO2 [Torr]", "pNO [Torr]", "T [K]"]
    data["rate N2 [mol/cm2·s"] = ammonia_oxidation([data[i] for i in data.columns[:3]], data["T [K]"])[0]
    data["rate NO [mol/cm2·s"] = ammonia_oxidation([data[i] for i in data.columns[:3]], data["T [K]"])[1]
    
    train,test = generate_excel(data, folder = r"Local\Data\Ammonia_Oxidation")
    return train,test
    
        
def propene_ammoxidation(C,T):
    """
    Reaction rate for the propene amoxidation in mol/m2·s.
    
    T is in Kelvin.
    
    Reference: https://www.sciencedirect.com/science/article/pii/S0021951716300239
    """
    pNH3, pC3H6, pO2, pACO, pH2O = C        # kPa
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
    
    rN2     = kN2*pNH3/(1 + kH2O*pH2O)
    rACN    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*(1/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5)) + kACON*pACO*pNH3
    rNBP    = 0.94*kRLS*pC3H6*(alpha*pNH3/(1 + alpha*pNH3))*((delta*pNH3 + epsilon*pO2**0.5)/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))
    rACO    = 0.94*kRLS*pC3H6*((1/(1 + alpha*pNH3))*(1/(1 + beta*pO2**0.5)) + (alpha*pNH3/(1 + alpha*pNH3))*(gamma*pH2O/(1 + gamma*pH2O + delta*pNH3 + epsilon*pO2**0.5))) - kACON*pACO*pNH3
    rOBP    = 0.94*kRLS*pC3H6*(1/(1 + alpha*pNH3))*((beta*pO2**0.5)/(1 + beta*pO2**0.5)) + 0.06*kRLS*pC3H6
    
    return rN2, rACN, rNBP, rACO, rOBP

def propene_ammoxidation_experiment(nsamples = 625):
    """
    We consider the initial total pressure in the paper. The temperature moves from 603 to 773 K.
    """
    pNH30   = 8.2           # kPa
    pC3H60  = 0.2           # kPa
    pO20    = 16.3          # kPa
    pACO0   = 0             # kPa
    pH2O0   = 2             # kPa
    lower = [pNH30*0.6, pC3H60*0.6, pO20*0.6, pACO0*0.6, pH2O0*0.6, 603]
    upper = [pNH30*2, (pC3H60+2)*2, pO20*2, (pC3H60+2)*2, pH2O0*4, 773]
    np.random.seed(1)
    lhs_data = lhs(6, nsamples)
    for column in range(lhs_data.shape[1]):
        lhs_data[:,column] = lhs_data[:,column]*(upper[column] - lower[column]) + lower[column]
    data = pd.DataFrame(lhs_data)
    data.columns = ["pNH3 [kPa]", "pC3H6 [kPa]", "pO2 [kPa]", "pACO [kPa]", "pH2O [kPa]", "T [K]"]
    data["rate N2 [mol/m2·s"]  = propene_ammoxidation([data[i] for i in data.columns[:5]], data["T [K]"])[0]
    data["rate ACN [mol/m2·s"] = propene_ammoxidation([data[i] for i in data.columns[:5]], data["T [K]"])[1]
    data["rate NBP [mol/m2·s"] = propene_ammoxidation([data[i] for i in data.columns[:5]], data["T [K]"])[2]
    data["rate ACO [mol/m2·s"] = propene_ammoxidation([data[i] for i in data.columns[:5]], data["T [K]"])[3]
    data["rate OBP [mol/m2·s"] = propene_ammoxidation([data[i] for i in data.columns[:5]], data["T [K]"])[4]

    train,test = generate_excel(data, folder = r"Local\Data\Propene_Ammoxidation")
    return train,test


if __name__ == "__main__":
    haber_bosch_experiment(625)
    ammonia_oxidation_experiment(625)
    propene_ammoxidation_experiment(625)