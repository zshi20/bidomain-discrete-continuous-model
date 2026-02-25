% This m-file creates time lapsed data in a single run. 

%% Read Files

current_folder = pwd;
addpath([current_folder, '\continuous domain']);
addpath([current_folder, '\discrete domain']);

[~, elem] = read_cells_ascii('mesh_data/circ_elements_finer2.txt');
[~, node] = read_nodes_ascii('mesh_data/circ_nodes_finer2.txt', ' ');

%%

elem = elem(:, 1 : 3);
node = node(:, [1 2]).*40; % now the domain is a unit circle of radius = 1 um

% scale
R = 1000;
node = node .* R; % circle of radius = 1000 um 

modelParam = struct;

% Temporal Parameters
modelParam.N_tstep = 288; % 24 hours
modelParam.tstep_size = 5; % min

% Perturbation Location
modelParam.perturb_mcrphge_loc = [0 0];

% Macrophage Death Params
modelParam.k_a = 0;  %min-1
modelParam.k_at = 0; %min-1 age-1

% Cell Transition Params
modelParam.k_cat = 1/120; %min-1
modelParam.K_d = 1; %[c]

% Monocyte adherance rate and steady state pop. density in mm2
modelParam.monocyte_density_SS = 20;           % mono count/mm2
modelParam.p_o = 1;                          % mono count/min/mm2
k_detach = modelParam.p_o/modelParam.monocyte_density_SS; %min^-1

% Run-N-Tumble Params
modelParam.monocyte_speed = 4; %um/min
modelParam.mcrphge_speed = 0.4; %um/min
modelParam.TaxisParam.l = 2e-5; %m
modelParam.TaxisParam.S = 5e6;

% Domain Params
modelParam.DomainParam.Dirichlet = freeBoundary(triangulation(elem,node));
modelParam.DomainParam.Neumann   = [];
modelParam.DomainParam.node      = node;
modelParam.DomainParam.elem      = elem;
modelParam.DomainParam.facet     = find_facets(elem);
modelParam.DomainParam.R         = R;

% Numerically Derived: Monocyte Leaving from edge Params
modelParam.MonocyteLeaveEdgeParam.mean   = 0.6; % numerically determined
modelParam.MonocyteLeaveEdgeParam.std    = 0.95; 

% Transport Params in the cases
% Case I: dispersed
% Case II: saturated
% Case III: clustered

%     CASE   I          II         III        
D_c_cases = [8.973E+05; 2.000E+04; 4.000E+03];
k_c_cases = [1.311E+00; 5.844E-01; 1.753E+00];
p_cases   = [7.914E+04; 6.000E+05; 2.048E+04];

case_ind = 3;
modelParam.sc_strength = p_cases(case_ind); % [c]*um2/min
modelParam.D_c =         D_c_cases(case_ind);
modelParam.k_c =         k_c_cases(case_ind);

%% Checking time scales

j01 = 2.4048;
TransportTimeConstant = exp(-((modelParam.D_c *(j01/R)^2 + modelParam.k_c)*modelParam.tstep_size));
TransportTimeConstantStd = 0.1;
IsSteadyTransport = TransportTimeConstantStd > TransportTimeConstant

%%
output = bidomain_model(modelParam, 'viz', 'off', 'temporalDataInConsole', 'on', 'temporalDataInOutput',  'on');

saveDir = 'bidomain temporal output';
save_simulation_data(saveDir, modelParam, output);

