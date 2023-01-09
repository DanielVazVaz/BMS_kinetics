# Streamlit app to perform these calculations
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import sys
BMS_FOLDER = r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\BMS"
sys.path.append(BMS_FOLDER)
sys.path.append(BMS_FOLDER + "\\" + "Prior")
from utils.kinetic_utils import PropeneAmmoxidation
import numpy as np
import os
from itertools import cycle

# Basic metadata
metadata = {
            "Local_path": r"C:\Users\dvazquez\Daniel\Articulos\Articulo_BMS_cineticas\Code\Local\Results",
            "Exp_name": r"PropAmox50points",
            "Initial_p": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
            "Exp_T": np.linspace(607, 773, 10),
            "W_values": np.linspace(0, 8, 50),
            "chemicals": ["C3H6", "N2", "ACN", "NBP", "ACO", "OBP", "NH3", "O2", "H2O"]
            }

test_conditions_dic = {#"T_experiments": np.linspace(613, 763, 10),
                       "W_values": np.linspace(0,8,50),
                       #"chemicals": ["ACN", "C3H6"],
                       #"integrator": "LSODA",
                       #"p0": np.array([8.0, 0, 0, 0, 0, 0, 8.0, 16.0, 10.5]),
                        }

## Configuration of the page
st.set_page_config(
    page_title="BMS kinetics",
    page_icon="🎈",
    layout = "wide"
)

# Expander with information
with st.expander("ℹ️ - About this app", expanded=True):

    st.write(
        """     
-   Only information for 50 points with and without noise, and 25 points with noise
-   Numerical instabilities are common
	    """
    )

    st.markdown("")
    
# SECTION 1: Choose the number of points, and if there is noise or not
st.markdown("# Choose the folder")
with st.container():
    col1, col2 = st.columns([1,1])
    with col1:
        folder_project = st.selectbox("Select folder", ["PropAmox25pointsNoisy", "PropAmox50points", "PropAmox50pointsNoisy"])
    with col2:
        noise    = st.selectbox("Noise in the measure", [True, False])
        
# SECTION 2: Choose the input conditions
st.markdown("# Choose the input conditions")
with st.form("my form"):
    with st.container():
        col1, col2, col3 = st.columns([1,1,1])
        with col1:
            bms_models = st.text_input("Excel for the BMS models", "BMS_results.xlsx")
        with col2:
            T_experiment = st.text_input("Test experiment temperature [K]", "650") 
        with col3:
            integrator   = st.selectbox("Select the integrator", ["LSODA", "RK45"])
        chemicals_chosen = st.multiselect("Choose the chemicals to use the BMS formula on", options = metadata["chemicals"], default = "ACN")
        chemical_chosen_plot = st.multiselect("Choose the chemicals to plot the BMS formula on", options = metadata["chemicals"], default = "ACN")
    with st.container():
        col = [0,0,0,0]
        initial_ps = [0]*9
        for row in range(2):
            col[row] = st.columns([1,1,1,1])
            for j in range(4):
                with col[row][j]:
                    initial_ps[4*row + j] = st.text_input(f'$p_{{{metadata["chemicals"][4*row +j]}}}^{{0}}$ [kPa]', str(metadata["Initial_p"][4*row + j]))
        initial_ps[8] = st.text_input(f'$p_{{{metadata["chemicals"][8]}}}^{{0}}$ [kPa]', str(metadata["Initial_p"][8]))
    initial_ps = [float(i) for i in initial_ps]
    initial_ps = np.array(initial_ps)
    radio = st.radio("What data to show in the table", ("Real", "True", "Calc"))
    submitted = st.form_submit_button("Submit")


# SECTION 3: Calculations and showing the table
st.markdown('# Table of values for the experiments')
# Reading data and performing alculations
if submitted:
    perform_exp = PropeneAmmoxidation(os.path.join(metadata["Local_path"], folder_project), noise=noise)
    perform_exp.calculate_test_experiment(BMS_models = bms_models, integrator = integrator, T_experiments=[float(T_experiment)],
                                          chemicals = chemicals_chosen, p0 = initial_ps,
                                        **test_conditions_dic)
    # Showing the table
    table = st.dataframe(perform_exp.test_experiments[radio])
    
# SECTION 4: Plot the results
st.markdown('# Plot')
if submitted:
    chemical_flows = ["F" + i + " [mol/s]" for i in metadata["chemicals"]]
    x_data = {}  # Will be a dictionary with key per specie, and the W values
    y_data = {}  # Will be a dictionary with key per specie, and the flows values
    SST = {}
    SSR = {}
    R2 = {}

    for chemical, chemical_flow in zip(metadata["chemicals"], chemical_flows):
        W_real = perform_exp.test_experiments["Real"].loc["Exp1", "W [m2]"]
        F_real = perform_exp.test_experiments["Real"].loc["Exp1", chemical_flow]
        
        F_true = perform_exp.test_experiments["True"].loc["Exp1", chemical_flow]
        notna_values = F_true.notna()
        F_true = F_true[notna_values]
        W_true = perform_exp.test_experiments["True"].loc["Exp1", "W [m2]"]
        W_true = W_true[notna_values]
        SST["True", chemical] = sum((F_true - F_real[notna_values])**2)
        SSR["True", chemical] = sum((F_real[notna_values] - F_real[notna_values].mean())**2)
        try:
            R2["True", chemical] = 1 - SST["True", chemical]/SSR["True", chemical]
        except:
            R2["True", chemical] = -300  
        
        
        F_calc = perform_exp.test_experiments["Calc"].loc["Exp1", chemical_flow]
        notna_values = F_calc.notna()
        F_calc = F_calc[notna_values]
        W_calc = perform_exp.test_experiments["Calc"].loc["Exp1", "W [m2]"]
        W_calc = W_calc[notna_values]
        SST["Calc", chemical] = sum((F_calc - F_real[notna_values])**2)
        SSR["Calc", chemical] = sum((F_real[notna_values] - F_real[notna_values].mean())**2)
        try:
            R2["Calc", chemical] = 1 - SST["Calc", chemical]/SSR["Calc", chemical]
        except:
            R2["Calc", chemical] = -300
        
        x_data["Real", chemical] = W_real 
        x_data["True", chemical] = W_true 
        x_data["Calc", chemical] = W_calc
        y_data["Real", chemical] = F_real 
        y_data["True", chemical] = F_true 
        y_data["Calc", chemical] = F_calc
        
    fig = go.Figure()
    color_cycler = cycle(['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52'])
    for chemical in chemical_chosen_plot:
        color = next(color_cycler)
        for state in ["Real", "True", "Calc"]:
            if state=="Real":
                label = chemical
                dash = None
                showlegend = True
                R2_label = ""
            elif state == "True":
                label = f"{chemical}_True"
                dash = "dot"
                showlegend = True
                R2_label =f"{R2['True', chemical]:.4f}"
            else: 
                label = f"{chemical}_Calc"
                dash = "dash"
                showlegend = True
                R2_label =f"{R2['Calc', chemical]:.4f}"
            fig.add_trace(go.Scatter(x = x_data[state, chemical], y = y_data[state, chemical], name = f"{label} ({R2_label})", legendgroup = chemical,
                                     showlegend = showlegend, mode = "lines", line = dict(color = color, dash = dash, width = 2)))

    
    fig.update_xaxes(
        #tickangle = 90,
        title_text = "W [m2]",
        title_font = {"size": 20},
        title_standoff = 25)

    fig.update_yaxes(
        title_text = "F [mol/s]",
        title_font = {"size": 20},
        title_standoff = 25)
    
    st.plotly_chart(fig, use_container_width= True)
        
