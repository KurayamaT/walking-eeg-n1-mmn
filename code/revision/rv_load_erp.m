function E = rv_load_erp(matpath)
% RV_LOAD_ERP  Load participant-average ERP arrays and return ROI indices.
%   Arrays are (time, channel, condition, participant); conditions are
%   sit, walk_free, walk_water, walk_clay.
S = load(matpath);
E.dev = S.dev_array; E.std = S.std_array; E.t = S.times_ms(:);
E.channels = cellstr(S.ch_names); E.conditions = cellstr(S.conditions);
E.subjects = cellstr(S.subjects);
assert(isequal(E.conditions(:)', {'sit', 'walk_free', 'walk_water', 'walk_clay'}));
assert(size(E.dev, 2) == 32 && size(E.dev, 3) == 4 && size(E.dev, 4) == 22);
roi = {'Fz', 'FC1', 'FC2', 'Cz'};
E.roi = cellfun(@(c) find(strcmp(E.channels, c), 1), roi);
E.n_std = S.n_std_mat; E.n_dev = S.n_dev_mat;
if size(E.n_std, 1) ~= 22, E.n_std = E.n_std'; E.n_dev = E.n_dev'; end
end
