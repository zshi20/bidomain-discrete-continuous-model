function vm_beta = find_cell_angle(l, K_d, S, grad_conc_X, conc_X)
%FIND_CELL_ANGLE finds the orientation of monocyte/macrophage from local
% IF concentration and concentration gradient. 
% Note: If kappa is smaller than 1e-5, the directionality of the distribution is
% negligible. But kappa --> 0 significantly increases the computational
% time, so we set the lower limit of kappa to 1e-5. 

vm_kappa = find_kappa(S, l, K_d, grad_conc_X, conc_X);

if vm_kappa < 1e-5
    vm_kappa = 1e-5;
end

vm_mu = atan2(grad_conc_X(2), grad_conc_X(1)); % angle corresponding to IF gradient
vm_beta = vmrand(vm_mu, vm_kappa);
end

