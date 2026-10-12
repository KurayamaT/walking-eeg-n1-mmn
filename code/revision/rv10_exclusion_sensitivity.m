function R = rv10_exclusion_sensitivity(erp)
% RV10_EXCLUSION_SENSITIVITY  N1 interaction after excluding the three participants
%   with the largest Sit-to-Walk-Free increase in late-window RMS, and leave-one-out
%   F values (Results 3.4).
E = rv_load_erp(erp);
dn = rv_roi_mean(E.dev, E, 80, 130); sn = rv_roi_mean(E.std, E, 80, 130);
y = (dn(:, 2) - dn(:, 1)) - (sn(:, 2) - sn(:, 1));
late = E.t >= 350 & E.t <= 450;
d = squeeze(mean(E.dev(:, E.roi, :, :), 2)); s = squeeze(mean(E.std(:, E.roi, :, :), 2));
q = (squeeze(sqrt(mean(d(late, :, :).^2, 1)))' + squeeze(sqrt(mean(s(late, :, :).^2, 1)))') / 2;
[~, ord] = sort(q(:, 2) - q(:, 1), 'descend'); keep = setdiff((1:22)', ord(1:3));
r = rv_one_sample(y(keep));
F = zeros(22, 1); p = zeros(22, 1);
for i = 1:22
    v = y; v(i) = []; [~, p(i), ~, st] = ttest(v, 0); F(i) = st.tstat^2;
end
R = table(string(strjoin(E.subjects(ord(1:3)), ' ')), r.n, r.F, r.p_t, min(F), max(F), min(p), max(p), sum(p < .05), ...
    'VariableNames', {'excluded', 'n', 'F', 'p', 'loo_F_min', 'loo_F_max', 'loo_p_min', 'loo_p_max', 'loo_significant'});
end
