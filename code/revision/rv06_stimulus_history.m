function R = rv06_stimulus_history(csvpath, qcpath)
% RV06_STIMULUS_HISTORY  Summaries of the stimulus-history check
%   (Supplementary Table S8; Supplementary Results S-R1).
%   qcpath: high-pass processing status CSV; defaults to recording_qc.csv
%   in the same folder as csvpath. The summary includes 1.0, 0.5 and 0.1 Hz,
%   with at least one retained trial and the recording-level channel exclusions.
%   csvpath: participant-level means by condition, filter and transition
T = readtable(csvpath, 'TextType', 'string');
if nargin < 2
    qcpath = fullfile(fileparts(csvpath), 'recording_qc.csv');
end
Q = readtable(qcpath, 'TextType', 'string');
excluded = Q(Q.n_bad_channels > 12 & ismember(Q.condition, ["sit", "walk_free"]) & ismember(Q.cutoff_hz, [1 .5 .1]), ...
    {'subject', 'condition', 'cutoff_hz', 'n_bad_channels'});
T.recording_eligible = true(height(T), 1);
for i = 1:height(excluded)
    m = T.subject == excluded.subject(i) & T.condition == excluded.condition(i) & abs(T.highpass_hz - excluded.cutoff_hz(i)) < 1e-9;
    T.recording_eligible(m) = false;
end
tr = ["S_to_S", "D_to_S", "S_to_D", "D_to_D"]; cond = ["sit", "walk_free"]; hp = [1 .5 .1];
S = table(); P = table();
for f = 1:numel(hp)
    for c = 1:2
        base = T.condition == cond(c) & abs(T.highpass_hz - hp(f)) < 1e-9;
        for g = 1:4
            rr = T(base & T.transition == tr(g), :); el = rr.recording_eligible & rr.retained_trials > 0;
            p1 = rv_mean_ci(rr.p1_mean_uv(el)); n1 = rv_mean_ci(rr.n1_mean_uv(el));
            S = [S; table(cond(c), hp(f), tr(g), min(rr.retained_trials(el)), max(rr.retained_trials(el)), ...
                min(rr.retained_trials), max(rr.retained_trials), p1.n, p1.mean, ...
                p1.ci_low, p1.ci_high, p1.positive, n1.n, n1.mean, n1.ci_low, n1.ci_high, 'VariableNames', ...
                {'condition', 'highpass_hz', 'transition', 'trials_min', 'trials_max', 'trials_min_all', 'trials_max_all', 'p1_n', 'p1_mean_uv', ...
                 'p1_ci_low_uv', 'p1_ci_high_uv', 'p1_positive', 'n1_n', 'n1_mean_uv', 'n1_ci_low_uv', 'n1_ci_high_uv'})]; %#ok<AGROW>
        end
        pr = {"standard", "S_to_S", "D_to_S"; "deviant", "S_to_D", "D_to_D"};
        for k = 1:2
            A = T(base & T.transition == pr{k, 2} & T.recording_eligible & T.retained_trials > 0, :);
            B = T(base & T.transition == pr{k, 3} & T.recording_eligible & T.retained_trials > 0, :);
            [~, ia, ib] = intersect(A.subject, B.subject);
            p1 = rv_mean_ci(B.p1_mean_uv(ib) - A.p1_mean_uv(ia)); n1 = rv_mean_ci(B.n1_mean_uv(ib) - A.n1_mean_uv(ia));
            P = [P; table(cond(c), hp(f), pr{k, 1}, pr{k, 3} + " minus " + pr{k, 2}, p1.n, p1.mean, p1.ci_low, ...
                p1.ci_high, n1.n, n1.mean, n1.ci_low, n1.ci_high, 'VariableNames', {'condition', 'highpass_hz', ...
                'current', 'contrast', 'p1_n', 'p1_mean_uv', 'p1_ci_low_uv', 'p1_ci_high_uv', 'n1_n', 'n1_mean_uv', ...
                'n1_ci_low_uv', 'n1_ci_high_uv'})]; %#ok<AGROW>
        end
    end
end
R.summary = S; R.paired = P;
R.eligibility = T; R.excluded_recordings = excluded;
end
