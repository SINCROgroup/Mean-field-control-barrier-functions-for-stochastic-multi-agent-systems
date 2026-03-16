function f = density_estimation(Z, X, dx, dy, a, N)

% Estimate a smooth periodic density from particle positions by combining
% a wrapped histogram with Gaussian convolution on the torus.

Nx = size(X, 1);

% Map particle coordinates to histogram bins in the periodic box.
ix = mod(floor((Z(:,1) + a) / (2*a) * Nx), Nx) + 1;
iy = mod(floor((Z(:,2) + a) / (2*a) * Nx), Nx) + 1;

% Build the empirical occupancy histogram.
rho = accumarray([iy, ix], 1, [Nx, Nx]);

% Convert counts into a density whose integral is approximately one.
rho = rho / (N * dx * dy);

% Build a smoothing kernel consistent with periodic convolution.
sigma = 0.5;
[Xg, Yg] = meshgrid((-Nx/2:Nx/2-1)*dx, (-Nx/2:Nx/2-1)*dy);
G = exp(-(Xg.^2 + Yg.^2)/(2*sigma^2));
G = fftshift(G);
G = G / sum(G(:));

% Smooth in Fourier space to preserve periodicity.
rho_smooth = ifft2(fft2(rho) .* fft2(G));
f = real(rho_smooth);

end
