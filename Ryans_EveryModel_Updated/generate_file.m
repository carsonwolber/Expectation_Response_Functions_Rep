addpath('DSGE_tools/');

%% RBC
%Load Parameters
[param,set] = rbc_parameters;

%Generate a mod object
[modl,param,set] = rbc_model(param,set);

%Generate model file
model_func(modl);

%% JR

%Load Parameters
[param,set] = jr_parameters;

%Generate a mod object
[modl,param,set] = jr_model(param,set);

%Generate model file
model_func(modl);

%% BS

%Load Parameters
[param,set] = bs_parameters;

%Generate a mod object
[modl,param,set] = bs_model(param,set);

%Generate model file
model_func(modl);

%% NK

%Load Parameters
[param,set] = nk_parameters;

%Generate a mod object
[modl,param,set] = nk_model(param,set);

%Generate model file
model_func(modl);

%% BLL

%Load Parameters
[param,set] = bll_parameters;

%Generate a mod object
[modl,param,set] = bll_model(param,set);

%Generate model file
model_func(modl);