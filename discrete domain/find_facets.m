function facets = find_facets(elem)
N_mesh_element = size(elem,1); % total # of mesh elements
facets = cell(N_mesh_element,1);
parfor i = 1:N_mesh_element
    facets{i} = nchoosek(elem(i,:),2); 
end
end