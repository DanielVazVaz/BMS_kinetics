# Repository for the BMS/Kinetics project
Here, we compile everything from the BMS/kinetics project so it is in a single place.

# Folder structure
The following folder structure is needed. The `main.py` script must be run from the PROJECT directory
```
PROJECT
|   README.md
|   .gitignore
|
+---scripts
|   |   main.py
|   |   interactive.ipynb
|
+---utils
|   |   utils.py
|   |   config.py
|   |   __init__.py
|
+---Local
|   +---BMS
|   |   |   All the BMS files
|   |
|   +---Data
|   |   +---ExperimentName01
|   |       |   data.xlsx  
|   |   +---ExperimentName02
|   |       |   data.xlsx   
|   |   +---ExperimentName03
|   |       |   data.xlsx  
|   |   +---... 
|   |
|   +---Results
```

This ensures that everything is imported correctly. Otherwise, you have to change import paths in both ```utils.py``` and ```main.py```. The ```config.py``` file indicates the path of the BMS files.

The position of the `Local/Data` and `Local/Results` folders is not that important, since these paths have to be provided to the main class to load data and save results, but it is a good practice to have them inside the `Local` folder. It also helps with the `.gitignore`.


# Structure of the `data.xls` files
The excel files require to have the data separated in training and test beforehand, using two diferent sheets called `train` and `test` respectively. The structure of both sheets is analogous. The first column must be an index column without any header, i.e., cell A1 must be empty, and cells A2,.., Ax, have to contain an index for the number of points (recommended from 0 to whichever amount of data you have). The rest of B1, C1, etc, have to contain the headers of the data, e.g., `Input1, Input2, Output1`. The cells below these headers contain the data points. 