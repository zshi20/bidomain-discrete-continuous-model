function [new_position, new_element] = cell_travel_in_a_dish(old_position, old_element, d, node, elem, facets)
% Given a complete description of the triangular mesh (with connectivity
% matrix "elem" and coordinate matrix "node"), the cell's current
% position and element, and the distance vector "d" which the cell is
% intended to travel toward, find out the new position and new element of
% the cell.

% Reference: Wang et al., An GPU-accelerated particle tracking method for 
% Eulerian–Lagrangian simulations using hardware ray tracing cores
% https://doi.org/10.1016/j.cpc.2021.108221

T_k = old_element;
X = old_position;
q = 1.000001;


while true

    % Initialize
    t_min = 1.1; 
    f_e = -1; 
    for f_i = 1:3 % check if the cell's intended trajectory (specified by vector d) would hit any facet of the element
        facet_xyz = node(facets{T_k}(f_i,:),:); 
        P_fi = facet_xyz(1,:);
        
        % Calculating the face normal (pointing outward)
        n_fi = [diff(facet_xyz(:,2)) -diff(facet_xyz(:,1))];
        apex_xyz = node(setxor(facets{T_k},facets{T_k}(f_i,:)),:);
        pv = (apex_xyz - P_fi);
        if dot(pv,n_fi) > 0
            n_fi = - n_fi;
        end

        t_fi = dot(P_fi-X, n_fi)/dot(d, n_fi);
        
%         if dot(d, n_fi) > 0 && (t_fi > -1e-4) && (t_fi < 1e-4) % edge case
%             t_min = t_fi;
%             f_e = f_i;
%             X_fe = X + q*t_fi*d; disp('edge case')
%         end

        if (t_fi < t_min) && (t_fi > 0) && (t_fi < 1) % find out which facet the cell will hit
            t_min = t_fi; 
            f_e = f_i; % update hit facet
            X_fe = X + q*t_fi*d; % where the particle hits the facet
        end
    end
    
    if f_e == -1 % particle stays in the cell
        X = X + d;
        break;
    else % particle hits facet
        X = X_fe;  
        T_k = find_neighbor_cell(T_k, elem, facets{T_k}(f_e,:));
        d = (1-t_min)*d;
    end

    if isempty(T_k) % if crosses a boundary
        break;
    end

end

new_position = X;
new_element = T_k;

end