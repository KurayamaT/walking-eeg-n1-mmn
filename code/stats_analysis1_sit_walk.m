function stats_analysis1_sit_walk(matpath)
% STATS_ANALYSIS1_SIT_WALK
%   Analysis 1 (Sit vs Walk-Free). Reproduces Table 2 from the deposited
%   subject-/condition-averaged ERP arrays. Frontocentral ROI {Fz, FC1, FC2, Cz};
%   N1 window 80-130 ms, MMN window 130-200 ms.
%
%   Usage: stats_analysis1_sit_walk('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
S = load(matpath);
dev = S.dev_array; sd = S.std_array; t = S.times_ms(:);
chs = cellstr(S.ch_names); conds = cellstr(S.conditions);

roi  = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx = cellfun(@(c) find(strcmp(chs, c), 1), roi);
n1w  = t >= 80  & t <= 130;
mmnw = t >= 130 & t <= 200;
csit = find(strcmp(conds, 'sit'), 1);
cw   = find(strcmp(conds, 'walk_free'), 1);

amp   = @(A, w) squeeze(mean(mean(A(w, ridx, :, :), 1), 2)).';   % -> (subj, cond)
devN1 = amp(dev, n1w);  stdN1 = amp(sd, n1w);
devM  = amp(dev, mmnw); stdM  = amp(sd, mmnw);
n = size(devN1, 1);

fprintf('\n=== Table 2A: N1-window cell means over ROI (mean +/- SD, uV), n=%d ===\n', n);
cm = @(v) sprintf('%+.3f +/- %.3f', mean(v), std(v));
fprintf('  sit-standard : %s   [target -0.140 +/- 0.401]\n', cm(stdN1(:, csit)));
fprintf('  sit-deviant  : %s   [target -0.764 +/- 0.736]\n', cm(devN1(:, csit)));
fprintf('  walk-standard: %s   [target -0.109 +/- 0.346]\n', cm(stdN1(:, cw)));
fprintf('  walk-deviant : %s   [target -0.453 +/- 0.376]\n', cm(devN1(:, cw)));

inter = (devN1(:, cw) - devN1(:, csit)) - (stdN1(:, cw) - stdN1(:, csit));
tI = mean(inter) / (std(inter) / sqrt(n)); F_I = tI^2; dzI = mean(inter) / std(inter);
pI = 1 - fcdf(F_I, 1, n - 1); etaI = F_I / (F_I + (n - 1)); pIw = signrank(inter);
fprintf('\n=== Table 2B: Stimulus x Condition interaction ===\n');
fprintf('  F(1,%d)=%.2f  p=%.3f  eta_p^2=%.3f  d_z=%.2f   [target F=6.26 p=.021 eta=.230 dz=0.53]\n', ...
        n - 1, F_I, pI, etaI, dzI);
fprintf('  Wilcoxon signed-rank on interaction contrast: p=%.3f   [target .028]\n', pIw);

devAvg = (devN1(:, csit) + devN1(:, cw)) / 2; stdAvg = (stdN1(:, csit) + stdN1(:, cw)) / 2;
dS = devAvg - stdAvg; tS = mean(dS) / (std(dS) / sqrt(n)); F_S = tS^2;
pS = 2 * (1 - tcdf(abs(tS), n - 1)); etaS = F_S / (F_S + (n - 1));
fprintf('\n=== Stimulus main effect (deviant vs standard) ===\n');
fprintf('  t(%d)=%+.2f  p=%.2g  F=%.2f  eta_p^2=%.3f   [target t=-5.76 p<.001 F=33.18 eta=.612]\n', ...
        n - 1, tS, pS, F_S, etaS);

dDev = devN1(:, cw) - devN1(:, csit); tD = mean(dDev) / (std(dDev) / sqrt(n));
pD = 2 * (1 - tcdf(abs(tD), n - 1)); dzD = mean(dDev) / std(dDev); pDw = signrank(dDev);
dStd = stdN1(:, cw) - stdN1(:, csit); tSt = mean(dStd) / (std(dStd) / sqrt(n));
pSt = 2 * (1 - tcdf(abs(tSt), n - 1));
fprintf('\n=== Table 2D: simple effects (Sit vs Walk-Free) ===\n');
fprintf('  deviant : t(%d)=%+.2f p=%.3f d_z=%.2f Wilcoxon p=%.3f   [target t=2.16 p=.042 dz=0.46 W=.046]\n', ...
        n - 1, tD, pD, dzD, pDw);
fprintf('  standard: t(%d)=%+.2f p=%.2f   [target t=0.46 p=.65]\n', n - 1, tSt, pSt);

mmnSit = devM(:, csit) - stdM(:, csit); mmnWk = devM(:, cw) - stdM(:, cw);
dM = mmnWk - mmnSit; tM = mean(dM) / (std(dM) / sqrt(n));
pM = 2 * (1 - tcdf(abs(tM), n - 1)); dzM = mean(dM) / std(dM);
fprintf('\n=== MMN amplitude Sit vs Walk-Free ===\n');
fprintf('  t(%d)=%+.2f p=%.3f d_z=%.2f   [target |t|=1.13 p=.270 dz=0.24]\n', n - 1, tM, pM, dzM);
end
