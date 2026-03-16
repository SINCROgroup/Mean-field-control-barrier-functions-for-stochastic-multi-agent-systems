function I = compute_interactions(Z, Zl, f1, f2, X, Y)

% For each follower, accumulate the contribution of all leaders through
% the precomputed periodic interaction kernels f1 and f2.
I = zeros(size(Z));
for i = 1:length(Z)
    rel_pos = angdiff(Zl, ones(size(Zl)).*Z(i, :));

    I1 = interp2(X, Y, f1, rel_pos(:, 1), rel_pos(:, 2));
    I2 = interp2(X, Y, f2, rel_pos(:, 1), rel_pos(:, 2));

    % Sum the pairwise contributions to obtain the net interaction drift.
    I(i, 1) = sum(I1);
    I(i, 2) = sum(I2);
end

end
