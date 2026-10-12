function q = rv_holm(p)
% RV_HOLM  Holm-Bonferroni adjusted p values.
p = p(:); m = numel(p);
[ps, order] = sort(p);
adj = cummax(min(1, (m - (1:m)' + 1) .* ps));
q = zeros(m, 1); q(order) = adj;
end
