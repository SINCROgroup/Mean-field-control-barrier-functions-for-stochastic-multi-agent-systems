function mean_obstacle = generate_obstacle_centers(N, L, radii, lambda)

% Randomly place the centers of the dangerous disks in the periodic box.
% The routine enforces a clearance from the goal region and a minimum
% pairwise separation between disks in the toroidal metric.

mean_obstacle = zeros(N, 2);
i = 0;

while i < N
    % Sample a candidate center uniformly in [-L,L]^2.
    x = -L + 2*L*rand;
    y = -L + 2*L*rand;

    % Keep the goal region around the origin free of dangerous disks.
    if norm([x y]) < 1 + radii(i+1)
        continue
    end

    % Check pairwise disk separation using periodic distances.
    ok = true;
    for j = 1:i
        dxy = wrapToL([x y] - mean_obstacle(j,:), L);
        d = norm(dxy);
        if d < lambda + radii(i+1) + radii(j)
            ok = false;
            break
        end
    end

    if ok
        % Accept the candidate center.
        i = i + 1;
        mean_obstacle(i,:) = [x y];
    end
end

end
