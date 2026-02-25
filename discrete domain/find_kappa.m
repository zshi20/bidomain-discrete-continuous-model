function kappa = find_kappa(S, l, K_d, gradient_c, c)
%FIND_KAPPA kappa = S*DFRO;
% S = Sensitivity, depending on cell type and chemoattractant
% DFRO = Difference in fractional receptor occupancy
% K_d = dissociation coefficient of ligand-receptor complex
% l = length of cell in the direction of chemotactic gradient

% Reference: Szatmary et al, 10.1016/j.jtbi.2017.05.014

dcdx = norm(gradient_c);
DFRO = l/K_d * dcdx * 1/(c/K_d + 1)^2; %um/M * M/um
kappa = S*DFRO;
end

