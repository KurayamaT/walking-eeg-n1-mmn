function T = rv03_covariate_slopes(erp)
% RV03_COVARIATE_SLOPES  Slopes of the N1 interaction contrast on three
%   mean-centred covariates (Supplementary Table S6; Results 3.4).
E = rv_load_erp(erp);
dn = rv_roi_mean(E.dev, E, 80, 130); sn = rv_roi_mean(E.std, E, 80, 130);
y = (dn(:, 2) - dn(:, 1)) - (sn(:, 2) - sn(:, 1));
late = E.t >= 350 & E.t <= 450;
d = squeeze(mean(E.dev(:, E.roi, :, :), 2)); s = squeeze(mean(E.std(:, E.roi, :, :), 2));
dr = squeeze(sqrt(mean(d(late, :, :).^2, 1)))'; sr = squeeze(sqrt(mean(s(late, :, :).^2, 1)))';
X = [(dr(:, 2) + sr(:, 2) - dr(:, 1) - sr(:, 1)) / 2, ...
     (dr(:, 2) - dr(:, 1)) - (sr(:, 2) - sr(:, 1)), ...
     (E.n_dev(:, 2) - E.n_dev(:, 1)) - (E.n_std(:, 2) - E.n_std(:, 1))];
lab = ["late_RMS_change"; "late_RMS_interaction"; "trial_count_interaction"];
T = table();
for j = 1:3
    A = [ones(22, 1), X(:, j) - mean(X(:, j))]; b = A \ y; res = y - A * b; df = 20;
    se = sqrt(diag(sum(res.^2) / df * inv(A' * A))); h = tinv(.975, df) * se; %#ok<MINV>
    rr = corrcoef(X(:, j), y);
    T = [T; table(lab(j), 22, b(2), b(2) - h(2), b(2) + h(2), 2 * tcdf(-abs(b(2) / se(2)), df), rr(1, 2), ...
        b(1), (b(1) / se(1))^2, 2 * tcdf(-abs(b(1) / se(1)), df), 'VariableNames', ...
        {'covariate', 'n', 'slope', 'ci_low', 'ci_high', 'slope_p', 'r', 'intercept', 'intercept_F', 'intercept_p'})]; %#ok<AGROW>
end
assert(all(abs(T.intercept - mean(y)) < 1e-12));
end
