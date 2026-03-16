function V = cbf_solver(N, nabla_delta_H, U, constraint_bound, alpha_val)

% Solve the quadratic program associated with the CBF filter:
% among all admissible velocities, choose the one closest to U.

vdim = 2*N;
u = U(:);
c = nabla_delta_H(:);

% Quadratic cost equivalent to minimizing ||v - u||^2.
H = speye(vdim);
f = -u;

% Encode the scalar CBF condition as a linear inequality for quadprog.
rhs = N * (constraint_bound - alpha_val);
A = -c.';
b = -rhs;

options = optimoptions('quadprog', 'Display', 'none');
v_star = quadprog(H, f, A, b, [], [], [], [], [], options);

% Reshape the optimizer output back into an N-by-2 velocity field.
V = reshape(v_star, N, 2);

end
