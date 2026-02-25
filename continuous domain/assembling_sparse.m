function A = assembling_sparse(node,elem,D,k)
%ASSEMBLING_SPARSE assembles the stiffness and mass matrix used for FEM
% computation (for the calculation of IF concentration). 
% This code is adapted from Chen's code (reference below). We
% restructured the code so that it is friendly to MATLAB's parallel
% computing.

% Reference: Chen L. Programming of finite element methods in matlab. 
% arXiv preprint arXiv:1804.05156. 2018 Apr 14.

N = size(node,1); % total number of nodes 
NT = size(elem,1); % number of elements
NV = size(elem,2); % number of vertices/nodes per mesh element
i = zeros(NT, NV*NV); j = zeros(NT, NV*NV); s = zeros(NT, NV*NV);

parfor t = 1:NT % for every element

    At = local_stiffness(node(elem(t,:),:),D) + local_reaction(node(elem(t,:),:),k);
    i_sliced = zeros(1, NV*NV);
    j_sliced = zeros(1, NV*NV);
    s_sliced = zeros(1, NV*NV);
    index_sliced = 0;

    for ti = 1:NV      % local/barycentric indices of a single elem

        for tj = 1:NV  % local/barycentric indices of a single elem

            index_sliced = index_sliced + 1; % sparce matrix data index
            i_sliced(index_sliced) = elem(t,ti); % i,j store the global indices of nodes, provided the local/barycentric indices
            j_sliced(index_sliced) = elem(t,tj);
            s_sliced(index_sliced) = At(ti,tj);

        end

    end

    i(t,:) = i_sliced;
    j(t,:) = j_sliced;
    s(t,:) = s_sliced;

end

A = sparse(i(:), j(:), s(:), N, N);       

end

