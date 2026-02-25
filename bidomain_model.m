
function output = bidomain_model(ModelParam, varargin)

% BIDOMAIN_MODEL simulates the post-perturbation response of inflammation
% via a discrete-continuum (bidomain) model.
%
%   Input: MODELPARAM is a struct that contains the model parameters. 
%
%          Optional: 'VIZ' ('on'/'off') is a binary parameter that specifies whether the immune
%          cell distribution and IF concentration colormap will be
%          displayed or not. Default setting is 'off', meaning it will not be
%          displayed. 
%
%          Optional: 'TEMPORALDATAINCONSOLE' ('on'/'off') is a binary parameter that specified
%          whether the population stats of the cells will be printed in
%          command line. Default setting is 'off', meaning it will not be printed.
%
%          Optional: 'TEMPORALDATAINOUTPUT' ('on'/'off') is a binary parameter that specified
%          whether the population stats of the cells in EVERY time step
%          will be outputted and subsequently saved. The default is 'off',
%          meaning ONLY the population stats of the last time step will be saved. 
%
%  Output: OUTPUT is a struct that contains the following arrays:
%
%          X is an N*2 array that contains the x, y coordinates of all cells at the last step. N
%          is the number of cells.
%
%          T is an N*1 array that contains the mesh element number of all cells at the last step.
%
%          cell_type is an N*1 array that records the identity of the cells at the last step. 
%          1 = monocyte, 2 = macrophage, 3 = cell not in domain. 
%
%          Optional: time_data is a struct that contains the above information in
%          every time step. This will only exist if 'TEMPORALDATAINOUTPUT'
%          is on. 

%% Parse input
p = inputParser;
addRequired(p, 'ModelParam', @isstruct);
addParameter(p, 'viz', 'off');
addParameter(p, 'temporalDataInConsole', 'off');
addParameter(p, 'temporalDataInOutput',  'off');

parse(p, ModelParam, varargin{:});


if ~strcmp(p.Results.viz, 'on') && ~strcmp(p.Results.viz, 'off')
    error('viz has to be either on or off')
elseif ~strcmp(p.Results.temporalDataInConsole, 'on') && ~strcmp(p.Results.temporalDataInConsole, 'off')
    error('temporalDataInConsole has to be either on or off')
elseif ~strcmp(p.Results.temporalDataInOutput, 'on') && ~strcmp(p.Results.temporalDataInOutput, 'off')
    error('temporalDataInOutput has to be either on or off')
else
    viz_on    = strcmp(p.Results.viz, 'on');
    disp_data = strcmp(p.Results.temporalDataInConsole, 'on');
    include_time_data = strcmp(p.Results.temporalDataInOutput, 'on');
end

%% Later used numbers/variables/functions go here
node    = ModelParam.DomainParam.node;
elem    = ModelParam.DomainParam.elem;
facets  = ModelParam.DomainParam.facet;

N_node          = size(node,1); % total # of nodes in the mesh
N_sd            = size(node,2); % # of dimensions
N_mesh_element  = size(elem,1); % total # of mesh elements

R = range(node(:,1),1)/2;
domain_area = pi*R^2; %um

% function handles
in_2D_triangle = @(X, T_k) inpolygon(X(1), X(2), node(elem(T_k,:),1), node(elem(T_k,:),2));

%% Simulation Parameters

% Boundary Conditions
Dirichlet   = ModelParam.DomainParam.Dirichlet; 
Newmann     = ModelParam.DomainParam.Neumann;

% MODEL PARAMETERS
% -----------------------------------------------------------
% Temporal Parameters
N_tstep                             = ModelParam.N_tstep;
tstep_size                          = ModelParam.tstep_size;

% Perturbation Location
assign_resident_mcrphge_location    = ModelParam.perturb_mcrphge_loc;

% Transport Params
D_diffusivity                       = ModelParam.D_c;
K_decay                             = ModelParam.k_c;
sc_strength                         = ModelParam.sc_strength; 

% Macrophage Death Params
k_a                                 = ModelParam.k_a;
k_at                                = ModelParam.k_at;

% Cell Transition Params
k_cat                               = ModelParam.k_cat;
K_d                                 = ModelParam.K_d; 

% Monocytes: Attachment, Steady state, and Detachment Params
monocyte_density_SS                 = ModelParam.monocyte_density_SS;
p_o                                 = ModelParam.p_o;
k_detach                            = p_o / monocyte_density_SS;

% Run-N-Tumble Params
monocyte_speed                      = ModelParam.monocyte_speed; 
mcrphge_speed                       = ModelParam.mcrphge_speed; 
l                                   = ModelParam.TaxisParam.l;
S                                   = ModelParam.TaxisParam.S;

% Monocytes leaving the edge Params
mean_monocyte_leave_edge_per_step   = ModelParam.MonocyteLeaveEdgeParam.mean; % numerically determined
std_monocyte_leave_edge_per_step    = ModelParam.MonocyteLeaveEdgeParam.std; 

% Number of cells
N_res_mcrphge            = size(assign_resident_mcrphge_location, 1);

N_monocyte_SS            = round(monocyte_density_SS * domain_area * 1e-6);

N_monocyte_per_step      = p_o * domain_area * tstep_size * 1e-6;
N_monocyte_per_step      = N_monocyte_per_step * (1-exp(-k_detach*tstep_size))/(k_detach*tstep_size);
N_monocyte_per_step      = round(N_monocyte_per_step);

N_edge_monocyte_per_step = round(normrnd(mean_monocyte_leave_edge_per_step, std_monocyte_leave_edge_per_step, [N_tstep, 1]));
N_edge_monocyte_per_step(N_edge_monocyte_per_step < 0) = 0;

% Total number of cells
N_random_cell = N_tstep*N_monocyte_per_step + N_monocyte_SS + N_res_mcrphge;
N_edge_cell   = sum(N_edge_monocyte_per_step);
N_cell = N_random_cell + N_edge_cell;

%% Initialize posit/concentrations and setting variables

T               = zeros(N_cell, 1);
X               = zeros(N_cell, N_sd);
C_cell          = zeros(N_cell, 1);
grad_C_cell     = zeros(N_cell, N_sd);
mcrphge_age     = zeros(N_cell, 1);

% Delete later
mono_ct = zeros(N_tstep,1);

T(1:N_random_cell, 1) = randi(N_mesh_element, [N_random_cell, 1]); 

for i = 1 : N_random_cell
    X(i,:) = random_X_in_elem(node,elem,T(i));
end

for i = N_random_cell + 1 : N_cell
    ang = rand*2*pi;
    X(i,:) = (0.95*R) .* [cos(ang), sin(ang)];
end

for i = N_random_cell + 1 : N_cell
    T(i,:) = find_element_from_xy_coordinates(X(i,:), node, elem);
end

t_counter       = 0;

cell_type       = zeros(N_cell,1);

time_data = struct;
time_data.X1 = zeros([N_cell, N_tstep]);
time_data.X2 = zeros([N_cell, N_tstep]);
time_data.cell_type = zeros([N_cell, N_tstep]);

%% Set cells

% Set steady state monocytes
cell_type  (1 : N_monocyte_SS) = 1;

zero_map = zeros(size(node,1),1);
if viz_on
    viz(node, zero_map, X(cell_type == 1, :), [], t_counter);
end

% spawn resident macrophages
cell_type (N_monocyte_SS + 1 : N_monocyte_SS + N_res_mcrphge) = 2;

if ~isempty(assign_resident_mcrphge_location)
    res_mcrphge_index = find(cell_type == 2);
    X(res_mcrphge_index, :) = assign_resident_mcrphge_location;

    for i = 1:N_res_mcrphge
        T(res_mcrphge_index(i)) = find_element_from_xy_coordinates(assign_resident_mcrphge_location(i,:), node, elem);
    end
end


%% Particle travels using Neighbor Searching Method


while true % loop thru every time steps
    % tic
    % ticBytes(gcp);
    
    % Identify secreting cells (sc)
    X_sc = X(cell_type==2,:);
    T_sc = T(cell_type==2,:);
    sc_magn = sc_strength.*ones(size(T_sc));

    % IF transport
    C_map = IF_transport(X_sc, T_sc, sc_magn, node, elem, [], Dirichlet, D_diffusivity, K_decay);

    if viz_on
        pause(0.1); 
        viz(node, C_map, X(cell_type==1,:), X(cell_type==2,:), t_counter);
    end

    % Identify in domain cells and create in domain variables (for more
    % efficient computing)
    in_dom_idx           = find(cell_type~=0);
    in_dom_grad_C_cell   = grad_C_cell(in_dom_idx, :);
    in_dom_C_cell        = C_cell(in_dom_idx, :);
    in_dom_X             = X(in_dom_idx, :);
    in_dom_T             = T(in_dom_idx, :);
    in_dom_cell_type     = cell_type(in_dom_idx, :);

    for j = 1:length(in_dom_idx)

        nodes_of_msh_element = unique(facets{in_dom_T(j,:)}); % node #s of the current hosting mesh element 
        
        % Find the concentration gradient of chemokines
        [grad_C_j, C_j] = find_grad_in_mesh_element(node(nodes_of_msh_element,:), ...        % node positions
                                                    in_dom_X(j,:), ...                                 % particle positions
                                                    C_map(nodes_of_msh_element));            % node concentrations
        in_dom_C_cell(j,:) = C_j;
        in_dom_grad_C_cell(j,:) = grad_C_j;

    end

    in_dom_X_next_step = in_dom_X;
    in_dom_T_next_step = in_dom_T;

    % Cell movement
    for j = 1:length(in_dom_idx) % loop thru every cell
        
        T_j = in_dom_T(j,:);
        X_j = in_dom_X(j,:);
        grad_C_j = in_dom_grad_C_cell(j, :);
        C_j = in_dom_C_cell(j, :);
        
        d0 = [];
        if in_dom_cell_type(j)==1
            d0 = monocyte_speed*tstep_size;
        elseif in_dom_cell_type(j)==2
            d0 = mcrphge_speed*tstep_size;
        else
            error('Cell does not have a class')
        end
           
        % Taxis (Von-Mises Distribution)
        beta_j = find_cell_angle(l, K_d, S, grad_C_j, C_j);
        d_j = d0.*[cos(beta_j), sin(beta_j)];
    
        % Cell move, provided the displacement vector
        [X_j, T_j] = cell_travel_in_a_dish(X_j, T_j, d_j, node, elem, facets);

        % Check if cells are in correct mesh elem
        if and(~in_2D_triangle(X_j, T_j),~isempty(T_j))
            % disp(['Particle position = ', num2str(X_k(1)), ', ', num2str(X_k(2))]);
            % error(['Particle is not in or on the edge of the element','-',num2str(T_k)]);
            T_j = find_element_from_xy_coordinates(X_j, node, elem);
            if T_j == -1
                in_dom_cell_type(j) = 0;
            end
        end

        % Record whether the cell is still in the domain or not
        if isempty(T_j)
            in_dom_cell_type(j) = 0;
            T_j = -2;
        end
        
        in_dom_X_next_step(j,:) = X_j;
        in_dom_T_next_step(j,:) = T_j;
    end
    
    in_dom_X = in_dom_X_next_step;
    in_dom_T = in_dom_T_next_step;

    % Update global variables
    grad_C_cell  (in_dom_idx,:)  = in_dom_grad_C_cell;
    C_cell       (in_dom_idx,:)  = in_dom_C_cell;
    X            (in_dom_idx,:)  = in_dom_X;
    T            (in_dom_idx,:)  = in_dom_T;
    cell_type    (in_dom_idx,:)  = in_dom_cell_type;
    
    % Monocyte maturation program
    is_mature = sample_from_probability((cell_type==1).*find_P_trans(k_cat, C_cell, tstep_size, K_d));
    cell_type(is_mature)          = 2;

    % Death program
    is_detached_monocyte   = sample_from_probability((cell_type==1).*(1-exp(-k_detach*tstep_size)));
    is_dead_mcrphge        = sample_from_probability(k_a*(cell_type==2) + k_at*mcrphge_age);
    
    cell_type(is_detached_monocyte | is_dead_mcrphge) = 0;

    % Monocytes enter the tissue
    monocyte_entry_idx = N_monocyte_per_step*t_counter + N_monocyte_SS + N_res_mcrphge + (1:N_monocyte_per_step)';
    if t_counter == 0
        monocyte_edge_entry_idx = N_random_cell + (1:N_edge_monocyte_per_step(1,:))';
    else
        monocyte_edge_entry_idx = N_random_cell + (sum(N_edge_monocyte_per_step(1:t_counter,:))+1:sum(N_edge_monocyte_per_step(1:t_counter+1,:)))';
    end
    cell_type([monocyte_entry_idx; monocyte_edge_entry_idx]) = 1;

    % Time counter ticking
    t_counter = t_counter + 1;
    mcrphge_age(cell_type==2) = mcrphge_age(cell_type==2) + tstep_size;
    
    current_time = tstep_size*t_counter; % min

    if disp_data
        disp(['Step ', num2str(t_counter), ...
            '; Time = ', num2str(floor(current_time/60)), ' hr ', ...
            num2str(mod(current_time, 60)), ' min', ...
          '; Monocyte count = ', num2str(sum(cell_type==1)), ...
          '; Macrophage count = ', num2str(sum(cell_type==2))]);
    end

    if include_time_data
        time_data.X1(:, t_counter) = X(:,1);
        time_data.X2(:, t_counter) = X(:,2);
        time_data.cell_type(:, t_counter) = cell_type;
        time_data.T(:, t_counter) = T;
    end

    monocyte_leave_ct = sum(T==-2);
    % mono_ct(t_counter) = sum(cell_type==1);

    % Check if there is any particle left in the domain
    if all(cell_type==0)
        disp('all cells have left the domain. Time stepping is terminated.')
        break;
    end

    % Stop the time stepping
    if t_counter > N_tstep-1
        % disp('Has reached the specified number of loops');
        break;
    end
    % tocBytes(gcp);
    % toc
   
end

if viz_on
    pause(0.1);
    viz(node, C_map, X(cell_type==1,:), X(cell_type==2,:), t_counter);
end

output = struct;
output.X = X;
output.T = T;
output.cell_type = cell_type;

if include_time_data
    output.time_data = time_data;
end
% output.monocyte_leave_ct = monocyte_leave_ct;
% output.mono_ct = mono_ct;
end

function dummy = viz(node, C_map, X_monocyte, X_macrophage, t_counter)

set(gcf,'position',[280,280,360,360])

visualize_C_map(node, C_map); hold on

if not(isempty(X_monocyte))
    s1 = scatter(X_monocyte(:,1), X_monocyte(:,2), 'yo', 'filled'); 
    s1.MarkerEdgeColor = 'k'; 
end

if not(isempty(X_macrophage))
    s2 = scatter(X_macrophage(:,1), X_macrophage(:,2),'ro', 'filled'); 
    s2.MarkerEdgeColor = 'k';

end

title(gca, ['Time Step = ', num2str(t_counter)], 'FontSize', 13, 'interpreter','latex');

hold off;
axis off;
end