function Fext = extendPeriodicField(F)
% Extend a 2-D periodic field sampled on a grid with duplicated endpoints.
%
% If x = -L:dx:L, the samples at -L and L represent the same point. The
% ghost sample at -L-dx must therefore copy L-dx, not the duplicated L.

    Fper = F;
    Fper(:, end) = Fper(:, 1);
    Fper(end, :) = Fper(1, :);

    Fext = [Fper(:, end-1), Fper, Fper(:, 2)];
    Fext = [Fext(end-1, :); Fext; Fext(2, :)];
end
