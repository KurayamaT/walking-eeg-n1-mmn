function check_topo_consistency(matpath)
% Verify the Table 3A "topographic consistency (Fisher-z LOO)" rows for Water vs Clay.
% Methods: N1 window uses the deviant-evoked and standard-evoked 32-ch maps separately;
% MMN window uses the deviant-minus-standard map. Topographic consistency = Fisher-z of
% each participant's map vs the leave-one-out group mean map (same condition).
if nargin < 1 || isempty(matpath); matpath = fullfile('..', 'data', 'erp_arrays_n22.mat'); end
S = load(matpath);
dev = S.dev_array; sd = S.std_array; t = S.times_ms(:);
conds = cellstr(S.conditions);
cw = find(strcmp(conds, 'walk_water'), 1); cc = find(strcmp(conds, 'walk_clay'), 1);
n1 = t >= 80 & t <= 130; mmn = t >= 130 & t <= 200;

% 32-ch maps (subj, ch) = mean over a window
mapof = @(A, w, c) squeeze(mean(A(w, :, c, :), 1)).';

fprintf('\n=== Table 3A topographic consistency (Fisher-z LOO), Water vs Clay ===\n');
report('N1 deviant-evoked map ', topoc(mapof(dev, n1, cw)),  topoc(mapof(dev, n1, cc)), '[target t=-1.07 p=.295]');
report('N1 standard-evoked map', topoc(mapof(sd,  n1, cw)),  topoc(mapof(sd,  n1, cc)), '[target t=+0.86 p=.402]');
report('MMN difference map    ', topoc(mapof(dev, mmn, cw) - mapof(sd, mmn, cw)), ...
                                  topoc(mapof(dev, mmn, cc) - mapof(sd, mmn, cc)), ...
                                  '[target Water 0.371 Clay 0.405 t=-0.28 p=.783 perm=.782 W=.949]');
end

function z = topoc(M)          % M: (subj, ch) -> Fisher-z LOO consistency per subject
n = size(M, 1); z = zeros(n, 1);
for s = 1:n
    loo = mean(M([1:s-1, s+1:n], :), 1);
    r = corr(M(s, :).', loo.');
    r = max(min(r, 0.999999), -0.999999);
    z(s) = atanh(r);
end
end

function report(name, zW, zC, tgt)
d = zW - zC; n = numel(d);
tt = mean(d) / (std(d) / sqrt(n)); pp = 2 * (1 - tcdf(abs(tt), n - 1)); dz = mean(d) / std(d);
perm = exact_signflip(d); wilc = signrank(d, 0, 'method', 'exact');
fprintf('  %s: Water %.3f+/-%.3f  Clay %.3f+/-%.3f  t(%d)=%+.2f p=%.3f d_z=%+.2f  perm=%.3f W=%.3f   %s\n', ...
        name, mean(zW), std(zW), mean(zC), std(zC), n - 1, tt, pp, dz, perm, wilc, tgt);
end


function p = exact_signflip(d)
% Exact two-sided sign-flip permutation p: enumerate all 2^n sign reassignments
% (deterministic, no RNG). Matches the exact permutation reported in Table 3A.
d = d(:); n = numel(d); N = 2^n; obs = abs(sum(d)); tol = 1e-9 * max(1, obs);
cnt = 0; batch = 2^16;
for st = 0:batch:N-1
    idx = (st:min(st + batch, N) - 1)';
    B = zeros(numel(idx), n);
    for b = 1:n; B(:, b) = bitget(idx, b); end
    cnt = cnt + sum(abs((1 - 2 * B) * d) >= obs - tol);
end
p = cnt / N;
end
