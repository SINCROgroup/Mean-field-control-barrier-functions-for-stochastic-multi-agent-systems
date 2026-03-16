close all
clear
clc

% MAIN SIMULATION SCRIPT FOR SHEPHERDING WITH CBF.
% It runs 50 Monte Carlo realizations, simulates the closed-loop dynamics,
% tracks safety/performance metrics, and stores trajectories for post-processing.

%% initializations
%Spatial meshes (to evaluate densities)
L = pi; %half-size of the domain in each direction
dx = 2*L/49;%0.05;%spatial step size
dy = 2*L/49;%0.05;
x = - L : dx : L; %x-mesh 1d
linspace(-L, L, 50);
y = - L : dy : L; %y-mesh 1d
[X, Y] = meshgrid(x, y); %2d mesh

%time mesh
dt = 0.01;
T = 100;
t = 0 : dt : T;

%parameters
D = 0.005;%diffusion coefficient
Lr = L;%length scale of repulsive interactions
%computation of interaction kernel
[q_T_x,q_T_y] = compute_q_T(X', Y',50, 2*L, 2*L, Lr);
f1 = -q_T_x';
f2 = -q_T_y';
f1 = [f1(:, end), f1, f1(:, 1)];
f1 = [f1(end, :); f1; f1(1, :)];
f2 = [f2(:, end), f2, f2(:, 1)];
f2 = [f2(end, :); f2; f2(1, :)];

%Desired follower density 
mu = 0;
nu = 0;
k1 = 9;
k2 = 9;
thetaX = (pi/L) * (X - mu);
thetaY = (pi/L) * (Y - nu);
rho_d = exp(k1*cos(X-mu)+k2*cos(Y-nu)+cos(X-mu)*cos(Y-nu)+sin(X-mu)*sin(Y-nu));
rho_d = rho_d ./ trapz(y, trapz(x, rho_d, 2)); %renormalize to followers number
rho_g = 3/sqrt(k1);

%Desired leader density 
load("desired_leaders_density2.mat")
rho_l_d = rho_bar_L;
A = 1 / (4 * L^2) * (0.12 - trapz(y, trapz(x, rho_l_d, 2)));
rho_l_d = rho_l_d + A;
rho_l_d_norm = rho_l_d ./ trapz(y, trapz(x, rho_l_d, 2));
Ml = 0.12;
Nl = 100;
Nf = 720;
N = Nl + Nf;
alpha = 1/N;

%meshes for periodic interpolation
xext = -L-dx : dx : L+dx;
yext = -L-dx : dy : L+dx;
[Xext, Yext] = meshgrid(xext, yext);

%Dangerous regions
C = 5;  %Number of disks
obstacles_radii = 0.5 * ones(C,1);
lambda = 0.25;   %Safe distance 
nk = 1000;   %Number of particles for each disk

%CBF parameters
constraint_bound = 0.001;
alpha_val = 0.1;
epsilon_Rk_L = 0.013;
epsilon_Rk_F = 0.01;

%Parameters for k(x,y)
sigma = 0.25;  
mu = 1;


%Path for saving the data 
current_path = pwd;
save_path = fullfile(current_path, 'el=0013-ef=001-5o/');

%Initial velocity of leaders and followers
sigma_micro = sqrt(D*dt*2);

%%Loop 
for kk = 1:50
    % Fix the seed so each realization is reproducible and identifiable.
    rng(kk);  

    % Sample dangerous disk centers and discretize the corresponding unsafe sets.
    mean_obstacle = generate_obstacle_centers(C, L, obstacles_radii, lambda);
    Xo = generate_obstacles(mean_obstacle, obstacles_radii, L, nk); 
    No = length(Xo);
    
    % Draw feasible initial conditions outside the dangerous-disk safety margins.
    Xl = generate_initial_positions(Nl, L, mean_obstacle, obstacles_radii, lambda);
    Xf = generate_initial_positions(Nf, L, mean_obstacle, obstacles_radii, lambda);
    
    % Estimate initial densities to initialize the KL-based control law.
    rho_l = density_estimation(Xl, X, dx, dy, L, Nl);
    rho_l_norm = rho_l ./ trapz(y, trapz(x, rho_l, 2));
    rho_f = density_estimation(Xf, X, dx, dy, L, Nf);
    rho_f_norm = rho_f ./ trapz(y, trapz(x, rho_f, 2)); 

    % Pointwise KL integrand used by the leader controller to match the target density.
    d_kl_norm = rho_l_norm .* log(rho_l_norm ./ rho_l_d_norm);

    % Preallocate trajectories and summary metrics for the whole simulation horizon.
    Xl_time = nan(length(t)*Nl, 2);
    Xf_time = nan(length(t)*Nf, 2);
    nf_inside = nan(size(t));
    nl_inside = nan(size(t));
    Hl = nan(size(t));
    Hf = nan(size(t));
    chi = nan(size(t));

    for i = 1 : length(t)

        % Store the current state before advancing one time step.
        Xl_time((i-1)*Nl+1:i*Nl, :) = Xl;
        Xf_time((i-1)*Nf+1:i*Nf, :) = Xf;

        if mod(i, 100) == 0
            fprintf('Iteration %d out %d\n', i, kk);
        end
        % Compute the nominal leader control from the macroscopic KL objective.
        U = leaders_control_macropt(Xl, x, y, L, Xext, Yext, Ml, d_kl_norm);%leader inputs
        % Compute the nominal follower drift induced by leader-follower interactions.
        Vlf = compute_interactions(Xf, Xl, 10*f1, 10*f2, Xext, Yext);%leader/follower interactions
        % The gain 10 accelerates convergence without changing the workflow.
        Vf_nom = alpha * Vlf;
 
        % Filter the leader input through the safety layer, then update positions.
        [Vl, Hl(i)] = compute_safe_velocity_Rk(Xl, Xo, Nl, U, sigma, mu, L, alpha_val, constraint_bound, epsilon_Rk_L, 0);
        Xl = wrapToL(Xl + dt * Vl, L);    
        rho_l = density_estimation(Xl, X, dx, dy, L, Nl);
        rho_l_norm = rho_l ./ trapz(y, trapz(x, rho_l, 2));

    
        % Filter the follower drift through the CBF and add Brownian noise.
        [Vf, Hf(i)] = compute_safe_velocity_Rk(Xf, Xo, Nf, Vf_nom, sigma, mu, L, alpha_val, constraint_bound, epsilon_Rk_F, D);
        Xf = wrapToL(Xf + dt * Vf + sigma_micro*randn(Nf,2), L);  % rumore sulla posizione
   
        % Refresh the KL integrand after updating the leader density.
        d_kl_norm = rho_l_norm .* log(rho_l_norm ./ rho_l_d_norm);
   
        % Compute minimum-image distances from every dangerous disk center.
        Df = sqrt((wrapToL(Xf(:,1) - mean_obstacle(:,1)', L)).^2 + (wrapToL(Xf(:,2) - mean_obstacle(:,2)', L)).^2);
        Dl = sqrt((wrapToL(Xl(:,1) - mean_obstacle(:,1)', L)).^2 + (wrapToL(Xl(:,2) - mean_obstacle(:,2)', L)).^2);

        % Count the agents currently lying inside at least one dangerous disk.
        nf_inside(i) = sum(any(Df < obstacles_radii', 2));
        nl_inside(i) = sum(any(Dl < obstacles_radii', 2));

        % Track how many followers have reached the goal region.
        Xf_norm = vecnorm(Xf, 2,2);
        chi(i) = sum(Xf_norm < rho_g);
    end
    
    % Save one file per realization for the plotting/figure scripts.
    folder_path = fullfile(save_path, num2str(kk));
    mkdir(folder_path);
    file_name = strcat('data', num2str(kk), '.mat');
    file_path = fullfile(folder_path, file_name);
    save(file_path, 'nl_inside', 'nf_inside', 'Hl', 'Hf', 'chi', 'Xl_time', 'Xf_time', '-mat')
end
