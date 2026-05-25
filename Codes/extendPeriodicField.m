function Fext = extendPeriodicField(F)
    Fext = [F(:, end), F, F(:, 1)];
    Fext = [Fext(end, :); Fext; Fext(1, :)];
end