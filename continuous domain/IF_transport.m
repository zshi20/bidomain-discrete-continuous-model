function u = IF_transport(X_sc, T_sc, sc_magn, node, elem, Neumann, Dirichlet, D, k)
% This code calcuates the IF concentration u using FEM.
% Reference: Chen L. Programming of finite element methods in matlab. 
% arXiv preprint arXiv:1804.05156. 2018 Apr 14.

N_node = size(node,1); % total # of nodes in the mesh

g_D = @(X) 0*ones(size(X,1),1);
g_N = @(X) 0*ones(size(X,1),1);

b = source_vector(X_sc, T_sc, sc_magn, node, elem);
A = assembling_sparse(node, elem, D, k);

%-------------------- Dirichlet boundary conditions------------------------
isBdNode = false(N_node,1);
isBdNode(Dirichlet) = true;
bdNode = find(isBdNode);
freeNode = find(~isBdNode);
u = zeros(N_node,1);
u(bdNode) = g_D(node(bdNode,:));
b = b - A*u;

%-------------------- Neumann boundary conditions -------------------------
if (~isempty(Neumann))
    Nve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
    edgeLength = sqrt(sum(Nve.^2,2));
    mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
    b = b + accumarray([Neumann(:),ones(2*size(Neumann,1),1)], ...
        repmat(edgeLength.*g_N(mid)/2,2,1),[N_node,1]);

end

% solve
u(freeNode) = A(freeNode,freeNode)\b(freeNode);

end

