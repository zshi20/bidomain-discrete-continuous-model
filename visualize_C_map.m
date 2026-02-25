function dummy = visualize_C_map(node, C_map)
%VISUALIZE2D Visualze the concentration map and cell position based on the
%inputted elem matrix (conn), node coordinate matrix (xyz_coord), 
%concentration that corresponds to each node (conc), and the cell position
%(X)

xmin = min(node(:,1)); xmax = max(node(:,1)); 
ymin = min(node(:,2)); ymax = max(node(:,2)); 
resolution = (xmax-xmin)/1000;
[xq,yq] = meshgrid(xmin:resolution:xmax, ymin:resolution:ymax);
cq = griddata(node(:,1),node(:,2),C_map,xq,yq);
h = pcolor(xq,yq,cq); shading flat; axis equal;
set(h, 'EdgeColor', 'none'); 
set(gca, 'ColorScale', 'linear');
cb = colorbar;
cb.TickDirection = 'in';
ylabel(cb,{' '},'FontSize',14,'Rotation',270);
clim([0 1])
colormap gray
end