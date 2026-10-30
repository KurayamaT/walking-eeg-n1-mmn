# Revision analyses (MATLAB)

MATLAB code for the analyses added during revision of the manuscript. Every function
works from participant-average ERP arrays or from derived tables; none reprocesses raw
EEG and none runs ICA. The canonical anchor is unchanged: N1 Sit × Walk-Free interaction
**F(1, 21) = 6.26**.

Requirements: MATLAB R2025b with the Statistics and Machine Learning Toolbox.

## Functions

| Function | Computes | Reported in |
|---|---|---|
| `rv01_filter_common_sample(paths)` | N1 interaction at three high-pass settings in the available sample and common 19 participants; paired changes between settings | Table 3 Panel B; Table S2 |
| `rv02_windows_multiplicity(erp)` | N1- and MMN-window tests for both contrasts, Holm–Bonferroni sensitivity, descriptive simple effects, mean interaction contrast with 95% CI | Table S3; Results 3.1; Fig. 2e |
| `rv03_covariate_slopes(erp)` | Slopes and 95% CIs of the N1 interaction contrast on three mean-centred covariates | Table S5; Results 3.4 |
| `rv04_latency_windows(erp)` | Additional latency-window sensitivity for both contrasts | Table S6; Results 3.2 |
| `rv05_early_positivity(erp, hpf_paths)` | 40–70-ms amplitudes, number of participants with positive means, leave-one-out range, peak latencies; early-waveform summaries in the common 19 participants | Table S7; Supplementary Results S-R1; Fig. S6 |
| `rv06_stimulus_history(csv, qc)` | Stimulus-history summaries and paired history contrasts at 1.0, 0.5 and 0.1 Hz, using the recording exclusions and retained-trial criterion in S-M10 | Table S8; S-R1 |
| `rv07_speed_order(csv)` | Self-selected walking speed and condition-by-position counts | Methods 2.2; Results |
| `rv08_stimulus_sequences(dir)` | Audit of the recorded standard/deviant sequences | Supplementary Methods S-M9 |
| `rv09_processing_counts(erp, csv, log)` | Condition-wise trial-rejection percentages, bad/interpolated channels, removed ICA components, boundary exclusions | Table S4 |
| `rv10_exclusion_sensitivity(erp)` | Interaction after excluding the three participants with the largest late-window RMS increase; leave-one-out F range | Results 3.4 |
| `rv11_filter_processing_status(csv)` | Processing outcomes of the high-pass reanalysis | Supplementary Methods S-M5 |

The `common18` output of `rv01` is retained for historical reproducibility and is not reported in the current manuscript.

Shared helpers: `rv_load_erp`, `rv_roi_mean`, `rv_one_sample`, `rv_shapiro_wilk`, `rv_holm`, `rv_mean_ci`.

Each function returns its results as MATLAB tables. Example:

```matlab
addpath('code/revision');
T = rv04_latency_windows('data/erp_arrays_n22.mat');
disp(T(:, {'component', 'lo_ms', 'hi_ms', 'contrast', 'estimate_uv', 'ci95_low_uv', 'ci95_high_uv', 'selected_p'}));
```

## Inputs

| Input | Used by |
|---|---|
| `erp_arrays_n22.mat` — canonical participant-average ERP arrays (in the data deposit) | `rv02`, `rv03`, `rv04`, `rv05`, `rv09`, `rv10` |
| ERP arrays reprocessed at 1.0, 0.5 and 0.1 Hz (three `.mat` files of the same layout) | `rv01`, `rv05` |
| Per-recording processing status of the high-pass reanalysis (`recording_qc.csv`) | `rv06`, `rv11` |
| Participant-level means by condition, filter and stimulus transition (`.csv`) | `rv06` |
| Event list per recording, `Sxx_<condition>_events.csv` (88 files) | `rv08` |
| Walking speed and condition order per participant (`.csv`) | `rv07` |
| Trigger and boundary counts per recording (`.csv`) and the preprocessing log (`.txt`) | `rv09` |

Only the canonical ERP arrays are part of the public data deposit. The other inputs are
derived, de-identified tables; they are not included in this repository.

### Stimulus-history inputs

Pass both paths explicitly when the files are in different folders:

```matlab
R = rv06_stimulus_history('data/stimulus_history_individual.csv', 'data/recording_qc.csv');
```

The participant-level CSV must include all three filter settings. The processing-status CSV
uses the columns `subject`, `condition`, `cutoff_hz` and `n_bad_channels`. If the second
argument is omitted, `recording_qc.csv` is read from the folder containing the first CSV.
Recordings with more than 12 bad channels are excluded at the corresponding setting;
each remaining transition requires at least one retained trial. For the reported data,
this gives n = 22 at 1.0 and 0.5 Hz, and n = 22 for Sit and n = 19 for Walk-Free at 0.1 Hz.
These are the existing eligibility rules in S-M10, not a new analysis.

## Conventions

- ROI: literature-defined frontocentral ROI {Fz, FC1, FC2, Cz}. Window boundaries are inclusive.
- Interaction contrast: (deviant Walk-Free − Sit) − (standard Walk-Free − Sit).
- Confidence intervals are Student-t intervals for arithmetic means.
- Test selection: one-sample t test when the Shapiro–Wilk p ≥ .05, otherwise the exact
  Wilcoxon signed-rank test. Shapiro–Wilk uses Royston's (1992) approximation.

See `VERIFICATION.md` for the archived numerical checks and the 30 October 2026 code-synchronisation record. Table and figure references above follow the current revised manuscript.
