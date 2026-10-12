function A = rv_roi_mean(X, E, lo, hi)
% RV_ROI_MEAN  Mean over the ROI electrodes and the inclusive window [lo, hi] ms.
%   Returns a participant-by-condition matrix.
w = E.t >= lo & E.t <= hi;
A = squeeze(mean(mean(X(w, E.roi, :, :), 1), 2)).';
end
