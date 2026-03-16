function X = generate_initial_positions(N, L, mean_obstacle, radii, epsilon)

% Rejection sampler for agent initialization in the periodic box.
% Candidate points are accepted only if they stay outside all dangerous
% disks enlarged by the safety margin epsilon.

thr = radii + epsilon; % Exclusion threshold for each dangerous disk

X = zeros(N,2);
ii = 0;

while ii < N
    % Sample uniformly in [-L,L]^2.
    pos = (2*L)*rand(1,2) - L;
    dxy = wrapToL(pos - mean_obstacle, L); % Periodic displacements
    dists = vecnorm(dxy,2,2); % Toroidal distances from dangerous disk centers

    if all(dists > thr)
        ii = ii + 1;
        X(ii,:) = pos;
    end
end

end
