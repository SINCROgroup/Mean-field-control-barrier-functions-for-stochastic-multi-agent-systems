function [rho_l_d, Mld] = leaders_density_reference(D, rho_d, x, y, f1, f2, L, dx, dy, N)
%given D, L, and rho_d, it computes the leaders density reference

[rho_d_x1, rho_d_x2] = gradient_2Dperiodic(rho_d, x, y, L);%feedforward action 

vfl_bar_x1 = N * D * rho_d_x1 ./ rho_d;
vfl_bar_x2 = N * D * rho_d_x2 ./ rho_d;

%rho_l_d = numerical_deconv(x,vfl_bar_x1, vfl_bar_x2, f1, f2);
rho_l_d = deconvolution_2D_pinv(f1, f2, vfl_bar_x1, vfl_bar_x2);

rho_l_d = gaussian_smoothing(rho_l_d,.5,x,y);
rho_l_d = rho_l_d - min([min(min(rho_l_d)),0]);
rho_l_d = circshift(circshift(rho_l_d,round(L/dx),2),round(L/dy),1);
%rho_l_d = rho_l_d(2:end-1,2:end-1);
Mld = trapz(x, trapz(y, rho_l_d, 2));
fprintf("Desired Leaders Mass = %2f \n",Mld);

end