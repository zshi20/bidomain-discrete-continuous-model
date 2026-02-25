function host_gradient = find_grad_of_node(host_coord, nghbr_coord, host_f, nghbr_f)
%FIND_GRAD_OF_NODE find the gradient (Ndim*1 matrix, with Ndim=2 or 3) of
% scalar function f at a hosting node.
% f(hosting_node) is known, and should be provided as an input to the fn. 
% Mechanism: 
% Solve for the least-sqr-solution x for the linear system b = Ax <==> 
% fi - f_host = (Xi - X_host) * gradient(f) at host, with i =
% 1,2,...N of nghbr_node

if ~isvector(nghbr_f) || ~isscalar(host_f)
    error('f(host) should be a scalar, and f(neighbor node) should be a vector, the length of which equals to the number of neighbors that the host has');
end

nghbr_f = reshape(nghbr_f, [], 1);
nghbr_n = length(nghbr_f);

if size(nghbr_coord,1) ~= nghbr_n
    error('incompatible sized matrices nghbr_coord and nghbr_f')
end

if size(host_coord,1) ~= 1
    error('host_coord is not correctly sized')
end

% b = Ax
b = nghbr_f - host_f;
A = nghbr_coord - repmat(host_coord,nghbr_n,1);
x = pinv(A)*b; % pinv = inv(A'*A)*A' or = A'*inv(A*A')

host_gradient = x';
end

