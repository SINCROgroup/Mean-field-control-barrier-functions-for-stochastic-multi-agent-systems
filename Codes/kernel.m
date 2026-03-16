function [K, dx, dy] = kernel(X, Y, L, sigma)
% Compute periodic pairwise displacements between agents X and sampled
% points Y belonging to the dangerous disks.
dx = wrapToL(X(:,1) - Y(:,1).', L);
dy = wrapToL(X(:,2) - Y(:,2).', L);

% Evaluate the isotropic Gaussian kernel on the torus.
D2 = dx.^2 + dy.^2;
K = exp(-D2 / (2*sigma^2));
end
