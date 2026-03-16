function [V, h_val] = compute_safe_velocity_Rk(X, Xo, N, U, sigma, mu, L, alpha, constraint_bound, epsilon_Rk, D)
% Compute a safe velocity field by filtering the nominal input U through a
% control barrier function based on the kernel overlap R_k.
%
% h = epsilon_Rk - R_k
% dh/dt + alpha(h) >= 0
%
% Inputs:
% X, Xo: agent positions and sampled points belonging to the dangerous disks
% N: number of agents
% U: nominal control input
% sigma, mu, L: kernel and scaling parameters
% alpha: class-K gain used in the barrier condition
% constraint_bound: extra safety margin in the QP constraint
% epsilon_Rk: admissible threshold for the overlap functional
% D: diffusion strength
%
% Outputs:
% V: safe velocity after CBF filtering
% h_val: current value of the barrier function

% Evaluate the Gaussian kernel between agents and the sampled dangerous disks.
[K, dx, dy] = kernel(X, Xo, L, sigma);

% Compute the kernel-weighted overlap with the dangerous disks.
Rk_val = mean(K(:));

% Gradient of the overlap functional with respect to the agent positions.
nabla_Rk = mu * [mean((dx ./ sigma^2) .* K, 2), ...
                 mean((dy ./ sigma^2) .* K, 2)];

% Barrier function h = epsilon_Rk - R_k.
h_val = epsilon_Rk - Rk_val;
alpha_val = alpha * h_val;

% Add the diffusion correction when the dynamics include noise.
if D > 0
    % Compute the Laplacian contribution of the kernel functional.
    r2 = dx.^2 + dy.^2;
    lap_mat = ((r2 - 2*sigma^2) / sigma^4) .* K;
    lap_Rk = -mean(lap_mat, 2);
    margin = D * mean(lap_Rk);
else
    margin = 0;
end

alpha_val = alpha_val + margin;

% Project the nominal velocity onto the CBF-admissible set.
V = cbf_solver(N, nabla_Rk, U, constraint_bound, alpha_val);
end
