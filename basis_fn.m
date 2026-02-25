function output = basis_fn(varargin)
% For triangular mesh element, linear basis function
% From "Construction of finite element spaces" pg. 110
% If two inputs to the function: they are 1. the triangle's x,y,z, which is
% a 3*3 matrix; 2. the query point's xyz in a 1*3 matrix
% If three inputs to the function: the additional input is the local node
% number which corresponds to the "A" in N_A. 

Ae = @(X1,X2,X3) 1/2*det([1 X1(1) X1(2); 1 X2(1) X2(2); 1 X3(1) X3(2)]); 

N1 = @(X1,X2,X3,X) 1/(2*Ae(X1,X2,X3)) * ...                     % 1/(2*Ae)
                 ((X2(1)*X3(2)-X3(1)*X2(2)) + ...               % (x2y3 - x3y2)
                 (X2(2)-X3(2))*X(1) + (X3(1)-X2(1))*X(2));      % (y2-y3)*x + (x3-x2)*y

X_verts = varargin{1};
X = varargin{2};
N_mat = [N1(X_verts(1,:), X_verts(2,:), X_verts(3,:), X);
         N1(X_verts(2,:), X_verts(1,:), X_verts(3,:), X);
         N1(X_verts(3,:), X_verts(2,:), X_verts(1,:), X)];

if length(varargin) == 2
    output = N_mat;
elseif length(varargin) == 3
    output = N_mat(varargin{3});
else
    error('Requiring two or three inputs')
end
end
