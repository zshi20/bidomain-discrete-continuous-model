function out_vector = source_vector(X_sc, sc_elementNo, sc_magn, node, elem)
% This code computed the forcing vector for FEM.

N_source = size(X_sc,1);
N_node = size(node, 1); %Number of nodes in the whole graph
n_node = size(elem, 2); %Number of nodes per element

if or(length(sc_elementNo)~=N_source, length(sc_magn)~=N_source)
    error('incompatibly sized input matrices')
end

out_vector = zeros(N_node,1);

% Note: assume col # in elem = barycentric ind
for k = 1:N_source
    for A = 1:n_node % Barycentric index, A = 1,2,3 for triangular mesh element
        X_vert = node(elem(sc_elementNo(k),:),:); % order is built in X_vert
        out_vector(elem(sc_elementNo(k),A)) = out_vector(elem(sc_elementNo(k),A)) + (sc_magn(k)*basis_fn(X_vert,X_sc(k,:),A));
    end
end

end