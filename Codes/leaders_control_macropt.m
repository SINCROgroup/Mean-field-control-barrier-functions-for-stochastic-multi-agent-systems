function [U, U_grid] = leaders_control_macropt(Xl, x, y, L, Xext, Yext, Ml, d_kl)

% Compute a periodic gradient-descent field for the KL objective that
% drives leaders toward the desired macroscopic density.
[u1, u2] = gradient_2Dperiodic(Ml*d_kl, x, y, L);

% Extend the field by one layer of cells so interpolation remains
% well-defined close to the periodic boundaries.
u1 = [u1(:, end), u1, u1(:, 1)];
u1 = [u1(end, :); u1; u1(1, :)];
u2 = [u2(:, end), u2, u2(:, 1)];
u2 = [u2(end, :); u2; u2(1, :)];

% Sample the feedback field at the current leader positions.
U1 = interp2(Xext, Yext, -100*u1, Xl(:, 1), Xl(:, 2));
U2 = interp2(Xext, Yext, -100*u2, Xl(:, 1), Xl(:, 2));

% Return both the grid field and the pointwise leader control.
U_grid = cat(3, -100*u1, -100*u2);
U = [U1, U2];

end
