function neighbor_cell_n = find_neighbor_cell(host_cell_n, connectivity, hit_facet)
% T_k = [] if find_neighbor_cell() can't find a neighbor cell that
% corresponds to the input given - which means the iteration was run in
% a border cell

cell_n = 1:size(connectivity,1);
is_neighborCell_or_hostCell = ones(length(cell_n),1); % initialze the bool array

for i = 1:length(hit_facet)
    is_neighborCell_or_hostCell = is_neighborCell_or_hostCell.*sum(connectivity == hit_facet(i),2);
end

is_hostCell = reshape(cell_n == host_cell_n, length(cell_n), []);
is_neighborCell = is_neighborCell_or_hostCell .* ~is_hostCell;
neighbor_cell_n = find(is_neighborCell);
end