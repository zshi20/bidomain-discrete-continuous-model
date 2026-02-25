% This file is an example of running 3D parameter sweep on the bidomain model.

%% Read Files

current_folder = pwd;
addpath([current_folder, '\continuous domain']);
addpath([current_folder, '\discrete domain']);

[~, elem] = read_cells_ascii('mesh_data/circ_elements_finer2.txt');
[~, node] = read_nodes_ascii('mesh_data/circ_nodes_finer2.txt', ' ');

%% Define constant variables

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
modelParam.p_o = 1;                            % mono count/min/mm2
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
modelParam.MonocyteLeaveEdgeParam.mean   = 0.6; 
modelParam.MonocyteLeaveEdgeParam.std    = 0.95; 

%% Define the parameter sweep

swpParam = struct;

swpParam.x.min = 1e3;
swpParam.x.max = 1e6;
swpParam.x.N   = 20;
swpParam.x.scale = 'log';
swpParam.x.paramName   = 'D_c';

swpParam.y.min = 0;
swpParam.y.max = 3;
swpParam.y.N   = 20;
swpParam.y.scale = 'linear';
swpParam.y.paramName  = 'k_c';

swpParam.z.min  = 10^3;
swpParam.z.max  = 10^(6.5);
swpParam.z.N    = 20;
swpParam.z.scale = 'log';
swpParam.z.paramName  = 'p';

swpParam.TransportTimeThreshold = 0.1;
fields = fieldnames(swpParam);

for i = 1:numel(fields)

    if ~isfield(swpParam.(fields{i}), 'min')
        continue;
    end

    if strcmp(swpParam.(fields{i}).scale, 'log')
        swpParam.(fields{i}).range = logspace(log10(swpParam.(fields{i}).min), log10(swpParam.(fields{i}).max), swpParam.(fields{i}).N);

    elseif strcmp(swpParam.(fields{i}).scale, 'linear')
        swpParam.(fields{i}).range = linspace((swpParam.(fields{i}).min), (swpParam.(fields{i}).max), swpParam.(fields{i}).N);
                                            
    else 
        error('scale has to be either linear or log (case sensitive)');

    end
end

[swpParam.x.grid, swpParam.y.grid, swpParam.z.grid] = ... 
           meshgrid(swpParam.x.range, swpParam.y.range, swpParam.z.range);

j01 = 2.4048;
TransportTimeConstant = exp(-((swpParam.x.grid *(j01/modelParam.DomainParam.R)^2 + swpParam.y.grid)*modelParam.tstep_size));
swpParam.IsSteadyTransport = swpParam.TransportTimeThreshold > TransportTimeConstant;

swpParam.x.steadyRange = swpParam.x.grid(swpParam.IsSteadyTransport);
swpParam.y.steadyRange = swpParam.y.grid(swpParam.IsSteadyTransport);
swpParam.z.steadyRange = swpParam.z.grid(swpParam.IsSteadyTransport);

%%
% Only simulate with the IF transport parameters that satisfies the steady
% state criterion
D_c_steadyRange = swpParam.x.steadyRange;
k_c_steadyRange = swpParam.y.steadyRange;
p_steadyRange   = swpParam.z.steadyRange;

outputMatrix = repmat(struct('X', [], 'T', [], 'cell_type', []), length(swpParam.x.steadyRange), 1);

parfor idx = 1:numel(D_c_steadyRange)

    tempModelParam = modelParam;
    tempModelParam.D_c           = D_c_steadyRange(idx);
    tempModelParam.k_c           = k_c_steadyRange(idx);
    tempModelParam.sc_strength   = p_steadyRange(idx);

    outputMatrix(idx) = bidomain_model(tempModelParam, 'viz', 'off');
end

% save the data
saveDir = 'bidomain output';
save_simulation_data(saveDir, modelParam, outputMatrix, swpParam);

