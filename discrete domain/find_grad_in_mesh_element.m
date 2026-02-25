function [grad_f_X0, f_X0] = find_grad_in_mesh_element(node_posit, particle_posit, f_node)

% Node_posit: N*M matrix; N = # of nodes, M = dimensions. 
%             The linear system is not going to solve unless N = M + 1
% Particle_posit: 1*M matrix; M = dimensions
% f_node: N*1 matrix

% Reference: https://math.stackexchange.com/questions/2627946/how-to-approximate-numerically-the-gradient-of-the-function-on-a-triangular-mesh/2632616#2632616

f_node = reshape(f_node, size(node_posit,1), []);
particle_posit = reshape(particle_posit, [], size(node_posit,2));

if size(f_node,2) ~= 1 || size(particle_posit,1) ~= 1
    error('Incompatible size of inputs, check input dimensions');
end

% Solve Ax = b

A = [node_posit' - repmat(particle_posit, size(node_posit,1), 1)'; ones(size(node_posit,1),1)']';

if ~(diff(size(A))==0)
    error('over or underconstrained linear system. A is not a square matrix')
end

b = f_node;
x = linsolve(A, b);
grad_f_X0 = x(1:end-1)';
f_X0 = x(end);
end