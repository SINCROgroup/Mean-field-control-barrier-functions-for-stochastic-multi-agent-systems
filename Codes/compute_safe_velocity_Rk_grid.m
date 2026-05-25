function [V, h_val] = compute_safe_velocity_Rk_grid(X, XX, YY, mu, constraint_bound, epsilon_Rk, D, U, Phi, Phi_x, Phi_y, N, lap_Phi, alpha)
% compute_safe_velocity_Rk_grid
% Computes a safe velocity field for N agents by enforcing a particle-level

    % Interpolate the gradient of Phi at the agent positions.
    grad_x = interp2(XX, YY, Phi_x, X(:,1), X(:,2));
    grad_y = interp2(XX, YY, Phi_y, X(:,1), X(:,2));

    % Gradient of R_k with respect to the particle positions.
    nabla_Rk = -mu * [grad_x, grad_y];

    % Evaluate the particle approximation of R_k by averaging Phi over agents.
    Phi_local = interp2(XX, YY, Phi, X(:,1), X(:,2));
    Rk_val = mean(Phi_local);

    % Barrier function: safety is encoded by h_val >= 0.
    h_val = epsilon_Rk - Rk_val;

    % Class-K term in the CBF condition.
    alpha_val = alpha * h_val;

    % Add the diffusion contribution, if present.
    if D > 0
        lap_local = interp2(XX, YY, lap_Phi, X(:,1), X(:,2));
        lap_Rk_val = -mean(lap_local);
        alpha_val = alpha_val + D * lap_Rk_val;
    end

    % Project the nominal velocity U onto the admissible half-space defined
    % by the CBF constraint. Since there is only one linear constraint, the
    % QP admits a closed-form projection and quadprog is not needed.
    u = U(:);
    c = nabla_Rk(:);
    rhs = N * (constraint_bound - alpha_val);

    cu = c.' * u;
    cn2 = c.' * c;

    % If the nominal input violates the CBF constraint, correct it along the
    % constraint normal direction. The threshold on cn2 avoids numerical
    % issues when the gradient is nearly zero.
    if cn2 > 1e-8 && cu < rhs
        u = u + ((rhs - cu) / cn2) * c;
    end

    % Reshape the corrected velocity vector back into an Nx2 matrix.
    V = reshape(u, N, 2);

end