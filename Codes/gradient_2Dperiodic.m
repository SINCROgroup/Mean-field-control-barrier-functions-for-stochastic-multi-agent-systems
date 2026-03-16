function [resx, resy] = gradient_2Dperiodic(A, x, y, L)

% Compute finite-difference gradients of a scalar field on a square
% periodic domain by tiling one additional period on every side.

x_vec = x;

% Extend the coordinates by one full period to the left and to the right.
x_vec_periodic = [x_vec - 2*L, x_vec, x_vec + 2*L];
y_vec_periodic = x_vec_periodic;

% Tile the field so the numerical gradient inherits periodic continuity.
A_periodic = [A, A, A; A, A, A; A, A, A];

[resx, resy] = gradient(A_periodic, x_vec_periodic, y_vec_periodic);

% Extract the gradient restricted to the central copy of the domain.
resx = resx(length(x_vec)+1 : 2*length(x_vec), length(x_vec)+1 : 2*length(x_vec));
resy = resy(length(x_vec)+1 : 2*length(x_vec), length(x_vec)+1 : 2*length(x_vec));

end
