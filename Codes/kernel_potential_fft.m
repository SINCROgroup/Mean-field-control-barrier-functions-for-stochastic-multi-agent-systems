function [Phi, Phi_x, Phi_y, lap_Phi] = kernel_potential_fft(rho_o_norm, X, Y, x, y, ~, sigma)
% Computes the convolution between the obstacle density rho_o_norm and a
% Gaussian kernel using FFTs. It also computes the corresponding spatial
% derivatives and Laplacian of the resulting potential.

    % Gaussian kernel centered at the origin.
    K_centered = exp(-(X.^2 + Y.^2)/(2*sigma^2));

    % Grid spacing and cell area for approximating the continuous integral.
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    cell_area = dx * dy;

    % Analytical derivatives of the Gaussian kernel.
    % These are consistent with the non-normalized kernel definition.
    Kx_centered = -(X / sigma^2) .* K_centered;
    Ky_centered = -(Y / sigma^2) .* K_centered;
    Klap_centered = ((X.^2 + Y.^2 - 2*sigma^2) / sigma^4) .* K_centered;

    % Move the kernel center to the first array entry, as required for
    % circular convolution using FFTs.
    K    = ifftshift(K_centered);
    Kx   = ifftshift(Kx_centered);
    Ky   = ifftshift(Ky_centered);
    Klap = ifftshift(Klap_centered);

    % Fourier transform of the obstacle density.
    F_rO = fft2(rho_o_norm);

    % Discrete approximation of the continuous convolution:
    % Phi(x) = int K(x-y) rho_o_norm(y) dy.
    % The multiplication by cell_area accounts for the quadrature weight.
    Phi      = real(ifft2(fft2(K)    .* F_rO)) * cell_area;
    Phi_x    = real(ifft2(fft2(Kx)   .* F_rO)) * cell_area;
    Phi_y    = real(ifft2(fft2(Ky)   .* F_rO)) * cell_area;
    lap_Phi  = real(ifft2(fft2(Klap) .* F_rO)) * cell_area;

end