function y = wrapToL(x, L)
% Map coordinates or displacements back into the periodic interval [-L, L).
y = mod(x + L, 2*L) - L;
end
