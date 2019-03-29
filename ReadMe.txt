=== README FOR CIRCADAPT MODEL AUGMENTATION INDEX SIMULATIONS ===

%% Running the model, generating the Data corresponding to data presented in the manuscript

To Run the CircAdapt Model for a single set of model parameters, execute and follow the instructions provided by the following script:

>>CircAdaptMain.m

To Run the Batch of CircAdapt simulations (corresponding to the parameter sweeps as detailed in the "Simulation Protocols" section), exectute the script:

>>CircAdaptBatchCall.m

On a commercial laptop, total simulation time is ~0.5 hrs, whereas at least 3.5 GB of free disk space is required
Please note that a folder with name "Simulations" will automatically be generated to contain the 153 (=17*9) simulations

%% Generating the plots and tables corresponding to those presented in the manuscript

To generate a figure or table as presented in the manuscript, evaluate the following scripts:

>>Generate_Fig4_Table1.m
>>Generate_Fig5ab.m
>>Generate_Fig5ab.m
>>Generate_Fig6.m
