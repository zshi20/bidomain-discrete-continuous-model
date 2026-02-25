function dummy = viz(node, C_map, X_monocyte, X_macrophage, t_counter)

set(gcf,'position',[280,280,360,360])

visualize_C_map(node, C_map); hold on

if not(isempty(X_monocyte))
    s1 = scatter(X_monocyte(:,1), X_monocyte(:,2), 'yo', 'filled'); 
    s1.MarkerEdgeColor = 'k'; 
end

if not(isempty(X_macrophage))
    s2 = scatter(X_macrophage(:,1), X_macrophage(:,2),'ro', 'filled'); 
    s2.MarkerEdgeColor = 'k';

end

title(gca, ['Time Step = ', num2str(t_counter)], 'FontSize', 13, 'interpreter','latex');

hold off;
axis off;
end