close all
clear
clc

% MAIN SCRIPT FOR THE COVERAGE APPLICATION.
% It simulates a population in a periodic domain with dangerous regions and stores
% trajectories together with dangerous regions incursion metrics for later plotting.



%time mesh
dt = 0.01;
T = 50;
t = 0 : dt : T;

%parameters
D = 0.05;%diffusion coefficient
Nf = 720;%Number of agents
L = pi; %half-size of the domain in each direction


%Dangerous disks
C = 5;  %Number of disks
obstacles_radii = 0.5*ones(C,1);
lambda = 0.25;   %Safe distance 
nk = 1000;   %Number of points for each disk
sigma_micro = sqrt(D*dt*2);

%CBF parameters
constraint_bound = 0.001;
alpha_val = 0.1;
epsilon_Rk_F = 0.01;

%Parameters for k(x,y)
sigma = 0.2;  
mu = 1;


%Path for saving the data 
current_path = pwd;
save_path = fullfile(current_path, 'heat-ef=001_5o');



%%Loop 
for kk = 1:50
    % Fix the seed so the Monte Carlo batch is reproducible.
    rng(kk);  

    % Sample dangerous disk centers and generate a dense point-cloud representation.
    mean_obstacle = generate_obstacle_centers(C, L, obstacles_radii, lambda);
    Xo = generate_obstacles(mean_obstacle, obstacles_radii, L, nk); 
    No = length(Xo);

    
    % Initialize agents uniformly in the unit disk around the origin.
    r = sqrt(rand(Nf,1));
    theta = 2*pi*rand(Nf,1);
    Xf = [r.*cos(theta) r.*sin(theta)];
    

    % Preallocate time histories and safety indicators.
    nf_inside = nan(size(t));
    Hf = nan(size(t));
    Vf_old = zeros(Nf, 2);
    a_max = 10;
    Xf_time = nan(length(t)*Nf, 2);
    for i = 1 : length(t)
        if mod(i, 100) == 0
            fprintf('Iteration %d out %d\n', i, kk);
        end

        % Store the current particle positions before propagation.
        Xf_time((i-1)*Nf+1:i*Nf, :) = Xf;
        
        % No nominal drift is prescribed in this experiment.
        U_nom = zeros(Nf,2);  
    
        % Apply the safety filter to the diffusive agents dynamics.
        [Vf, Hf(i)] = compute_safe_velocity_Rk(Xf, Xo, Nf, U_nom, sigma, mu, L, alpha_val, constraint_bound, epsilon_Rk_F, D);
        
        Xf = wrapToL(Xf + dt * Vf + sigma_micro*randn(Nf,2), L);  % rumore sulla posizion
        Vf_old = Vf;
    
        % Compute minimum-image distances from every dangerous disk center.
        Df = sqrt((wrapToL(Xf(:,1) - mean_obstacle(:,1)', L)).^2 + (wrapToL(Xf(:,2) - mean_obstacle(:,2)', L)).^2);
    
        % Count agents inside at least one dangerous-disk region.
        nf_inside(i) = sum(any(Df < obstacles_radii', 2));
    end
    
    % Save one file per realization for the figure-generation scripts.
    folder_path = fullfile(save_path, num2str(kk));
    mkdir(folder_path);
    file_name = strcat('data', num2str(kk), '.mat');
    file_path = fullfile(folder_path, file_name);
    save(file_path, 'nf_inside', 'Hf', 'Xf_time', '-mat')
end
