function R = rv11_filter_processing_status(qc_csv)
% RV11_FILTER_PROCESSING_STATUS  Processing outcomes of the high-pass reanalysis
%   (Supplementary Methods S-M5; Results 3.2). qc_csv: one row per recording and
%   cutoff with the processing status, the number of flagged channels and whether
%   FastICA met its convergence criterion.
T = readtable(qc_csv, 'TextType', 'string');
cut = [1 .5 .1]; S = table();
for f = 1:3
    r = T(abs(T.cutoff_hz - cut(f)) < 1e-9, :); bad = r.status ~= "ok";
    S = [S; table(cut(f), height(r), sum(~bad), sum(bad & startsWith(r.error, "too many bad ch")), ...
        sum(bad & ~startsWith(r.error, "too many bad ch")), 'VariableNames', ...
        {'cutoff_hz', 'recordings', 'erp_generated', 'bad_channel_limit', 'other_failure'})]; %#ok<AGROW>
end
m = T(abs(T.cutoff_hz - .1) < 1e-9 & T.condition == "walk_free" & T.status ~= "ok", {'subject', 'n_bad_channels', 'error'});
i = T(T.status == "ok" & lower(string(T.ica_converged)) == "false", {'cutoff_hz', 'subject', 'condition'});
R.status = S; R.missing_walk_free_0p1 = sortrows(m, 'subject'); R.ica_iteration_limit = sortrows(i, {'cutoff_hz', 'subject'}, {'descend', 'ascend'});
end
