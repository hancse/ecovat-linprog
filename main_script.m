%% Prepare the workspace

clear all
close all
clc

%% Load ECOVAT buffer physical parameters:
ECOVAT_Parameters;

%% Load the external datasets used in the optimization routine:
External_datasets;

%% Create the constriants, objective functions and solve:
Problem_formulation_oct;

%% Post-process and data visualization: