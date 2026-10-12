function R = rv09_processing_counts(erp, denominators_csv, logpath)
% RV09_PROCESSING_COUNTS  Condition-wise processing counts of the primary analysis
%   (Supplementary Table S5). Rejection percentage = amplitude-rejected epochs
%   among candidate epochs lying fully within the recording.
%   erp: canonical arrays (retained trial counts); denominators_csv: trigger and
%   boundary counts per recording; logpath: preprocessing log (bad channels, removed ICs).
E = rv_load_erp(erp);
D = readtable(denominators_csv, 'TextType', 'string');
txt = string(splitlines(fileread(logpath)));
cond = ["sit", "walk_free", "walk_water", "walk_clay"]; subj = string(E.subjects(:));
P = table();
for i = 1:height(D)
    si = find(subj == D.subject(i)); ci = find(cond == D.condition(i));
    tok = regexp(txt, "^\s*OK\s+" + D.subject(i) + "\s+" + D.condition(i) + "\s+n_std=\s*(\d+)\s+n_dev=\s*(\d+)\s+bad_ch=(.*?)\s+ICs_rm=(\d+)", 'tokens', 'once');
    tok = tok(~cellfun(@isempty, tok)); assert(numel(tok) == 1); tok = tok{1};
    assert(str2double(tok(1)) == E.n_std(si, ci) && str2double(tok(2)) == E.n_dev(si, ci));
    b = strtrim(tok(3));
    if startsWith(b, "zeros"), nbad = 0; else, nbad = numel(str2num(b)); end %#ok<ST2NM>
    for st = ["standard", "deviant"]
        total = D.(st + "_total")(i); bnd = D.(st + "_boundary")(i); valid = total - bnd;
        if st == "standard", kept = E.n_std(si, ci); else, kept = E.n_dev(si, ci); end
        assert(kept >= 0 && kept <= valid);
        P = [P; table(D.subject(i), D.condition(i), st, total, bnd, valid, kept, valid - kept, ...
            100 * (valid - kept) / valid, nbad, str2double(tok(4)), 'VariableNames', {'subject', 'condition', ...
            'stimulus', 'total_triggers', 'boundary_excluded', 'boundary_valid', 'retained', 'amplitude_excluded', ...
            'amplitude_rejection_pct', 'bad_interpolated_channels', 'removed_ICs'})]; %#ok<AGROW>
    end
end
S = table();
for c = cond
    for st = ["standard", "deviant"]
        r = P(P.condition == c & P.stimulus == st, :);
        S = [S; table(c, st, height(r), sum(r.total_triggers), sum(r.boundary_excluded), sum(r.retained), ...
            mean(r.amplitude_rejection_pct), std(r.amplitude_rejection_pct), mean(r.bad_interpolated_channels), ...
            std(r.bad_interpolated_channels), mean(r.removed_ICs), std(r.removed_ICs), 'VariableNames', ...
            {'condition', 'stimulus', 'n', 'total_triggers', 'boundary_excluded', 'retained', 'mean_rejection_pct', ...
             'sd_rejection_pct', 'mean_bad_interpolated', 'sd_bad_interpolated', 'mean_removed_ICs', 'sd_removed_ICs'})]; %#ok<AGROW>
    end
end
R.participants = P; R.summary = S;
end
