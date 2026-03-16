function Xo = generate_obstacles(mean_obstacle, r, L, nk)

% Sample nk points uniformly inside each dangerous disk.
% This point cloud is the discrete representation of the dangerous set
% used by the kernel-based barrier function.

M = size(mean_obstacle,1);
Xo = zeros(nk*M, 2);
ofs = 0;

for k = 1:M
    th = 2*pi*rand(nk,1);
    rr = r(k)*sqrt(rand(nk,1)); % Uniform density in a disk
    pts = [mean_obstacle(k,1) + rr.*cos(th), mean_obstacle(k,2) + rr.*sin(th)];

    % Wrap points back into the periodic domain if a disk crosses the boundary.
    pts = mod(pts + L, 2*L) - L;

    Xo(ofs+(1:nk), :) = pts;
    ofs = ofs + nk;
end

end
