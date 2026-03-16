function [q_T_x,q_T_y] = compute_q_T(grid_x,grid_y,N_sum_q,L_x,L_y,Lr)
%% This function computes the discretized version of the interaction kernel.
% The kernel is evaluated through a periodic image summation, so it is
% consistent with the toroidal domain used in the simulations.
%
% Input:
% grid_x: Mesh grid (x coordinate) of the cell centers
% grid_y: Mesh grid (y coordinate) of the cell centers
% N_sum_q: Number of periodic images used in each direction
% L_x: Domain length along the x-axis
% L_y: Domain length along the y-axis
% Lr: Interaction length scale
%
% Output:
% q_T_x: Discretized kernel in the x-direction
% q_T_y: Discretized kernel in the y-direction

N_xC = size(grid_x,1); % Number of cells along the x-axis
N_yC = size(grid_x,2); % Number of cells along the y-axis

q_T_x = zeros(N_xC,N_yC); % Preallocation
q_T_y = zeros(N_xC,N_yC); % Preallocation

%% Computation
i = round(N_xC/2);
j = round(N_yC/2);
x_i = grid_x(i,1);
y_j = grid_y(1,j);

% Evaluate the interaction of each grid point with the reference cell at
% the center of the computational domain.
for k = 1:N_xC
    x_k = grid_x(k,1);
    for l = 1:N_yC
        y_l = grid_y(1,l);
        for m = -N_sum_q:N_sum_q
            for n = -N_sum_q:N_sum_q
                if k==i && l==j && m==0 && n==0
                    % Skip the self-interaction term because it is singular.
                else
                    theta = atan2(y_j-y_l+L_y*n, x_i-x_k+L_x*m);
                    dist = sqrt((x_i-x_k+L_x*m)^2 + (y_j-y_l+L_y*n)^2);
                    q_T_x(k,l) = q_T_x(k,l) + exp(-dist/Lr)*cos(theta);
                    q_T_y(k,l) = q_T_y(k,l) + exp(-dist/Lr)*sin(theta);
                end
            end
        end
    end
end
