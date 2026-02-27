function which_elem = find_element_from_xy_coordinates(coords, node, elem)
%FIND_ELEMENT_FROM_XY_COORDINATES this function finds the mesh element number that
%corresponds to the provided x, y coordinates.
%coords:     1x2 array of (x, y) coordinates (cant be more than one set of coordinate!)
%node:       Nx2 node coordinate matrix
%elem:       Mx3 connectivity matrix
%which_elem: The index of the element the coordinate belongs into

% A temporary variable (3x2 matrix) that contains the vertice coordinate
% of the current triangle in the loop
    triangle_node = zeros(3,2);

    for i = 1:size(elem, 1) % Going through all mesh elements

        % Get the coordinates of the vertices of the current triangle element
        for j = 1:3
            triangle_node(j,:) = node(elem(i,j),:);
        end

        % Calculate Jacobian
        J = [triangle_node(2,:)-triangle_node(3,:); 
             triangle_node(1,:)-triangle_node(3,:)]';

        bary_node = J\(coords-triangle_node(3,:))'; % Calculate the barycentric corrdinate

        if bary_node(1) > 0 && bary_node(2) > 0 && bary_node(1) + bary_node(2) < 1  % if within the unit triangle (bary coord)
            which_elem = i;
            return;
        end      
    end

    % If no element is found, return -1
    which_elem = -1;
end




