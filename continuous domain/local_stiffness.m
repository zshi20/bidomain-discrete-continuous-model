

function out_matrix = local_stiffness(varargin)

if length(varargin) == 2
    X_vert = varargin{1};
    D = varargin{2};
elseif length(varargin) == 1
    X_vert = varargin{1};
    D = 1;
else
    error('There cannot be more than 2 inputs or less than 1 input')
end

nsd = size(X_vert,2);
% Check D dimension
if isscalar(D)
    D = eye(nsd).*D; % I*D
elseif ~isequal(size(D),[nsd,nsd])
    error('Diffusion coefficient must be a square matrix or a scalar (assuming isotropic)')
end


% X_vert <- 3*2
% X_vert = [x_A1, y_A1;
%           x_A2, y_A2;
%           x_A3, y_A3]

J = [X_vert(1,:)-X_vert(3,:); X_vert(2,:)-X_vert(3,:)]';

gradN =   [1   0;  % d/dksi(N_A1)
           0   1;  % d/dksi(N_A2)
          -1  -1]; % d/dksi(N_A3)

integral_factor = 1/2; 
%out_matrix = gradN*pinv(J)*transpose(gradN*pinv(J))*norm(cross(J(:,1),J(:,2)))*integral_factor;
if size(J,1) == size(J,2)
    out_matrix = gradN/J*D*transpose(inv(J))*transpose(gradN)*abs(det(J))*integral_factor;
else
    out_matrix = gradN*pinv(J)*D*transpose(gradN*pinv(J))*norm(cross(J(:,1),J(:,2)))*integral_factor;
end

end