

function out_matrix = local_reaction(varargin)

if length(varargin) == 2
    X_vert = varargin{1};
    k = varargin{2};
elseif length(varargin) == 1
    X_vert = varargin{1};
    k = 1;
else
    error('There cannot be more than 2 inputs or less than 1 input')
end

% X_vert <- 3*2
% X_vert = [x_A1, y_A1;
%           x_A2, y_A2;
%           x_A3, y_A3]

J = [X_vert(1,:)-X_vert(3,:); X_vert(2,:)-X_vert(3,:)]';

out_matrix = 1/24 .* [2   1   1;  
                      1   2   1;  
                      1   1   2]; 

if size(J,1) == size(J,2)
    out_matrix = out_matrix.*(k*abs(det(J)));
else
    out_matrix = out_matrix.*(k*norm(cross(J(:,1),J(:,2)))); 
end
end