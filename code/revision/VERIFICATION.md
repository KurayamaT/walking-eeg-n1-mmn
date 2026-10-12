# Verification of the revision analyses (12 October 2026)

The functions in this folder were run on the authors' analysis machine and compared with
(1) the outputs archived when the revision analyses were first computed and (2) the values
reported in the revised manuscript, Supporting Information and response letter.

Two kinds of agreement are distinguished below.

- **Agreement within a numerical tolerance**: unrounded MATLAB values against the unrounded
  archived values; the absolute difference is at or below the stated tolerance.
- **Agreement at the reported precision**: MATLAB values formatted to the number of decimals
  and the notation used in the manuscript, compared as strings with the reported cells.

## Summary

| Check | Kind of agreement | Size | Result |
|---|---|---|---|
| Archived outputs of the revision analyses | numerical tolerance | 23 comparisons, 3,587 values | all agree |
| Reported values (tables and text) | reported precision | 374 cells | all agree |
| The ten scripts in `code/`, rerun | reported precision (output read against the targets printed by the scripts and the manuscript) | 10 scripts | agree |

No numerical discrepancy was found.

## What was computed

Eleven groups of analyses (`rv01`–`rv11`): high-pass filter comparison in common samples;
processing status of the high-pass reanalysis; N1- and MMN-window tests with Holm–Bonferroni
sensitivity; simple effects and the interaction contrast; covariate slopes; additional latency
windows; early positivity; stimulus history; walking speed and condition order; stimulus
sequences; processing counts; exclusion sensitivity. The README lists where each is reported.

## Environment

- MATLAB R2025b Update 5 (25.2.0.3177638), maca64
- Statistics and Machine Learning Toolbox 25.2 (R2025b)
- macOS 26.5.2 (arm64)
- Run on 12 October 2026

## Input files

| Input | Bytes | SHA-256 |
|---|---|---|
| `erp_arrays_n22.mat` — canonical participant-average ERP arrays (in the public data deposit) | 33,438,472 | `3e8bbe957367fb85cc91f1c888f65201fd48af8ec4bf64504bc7363e0eb58beb` |
| ERP arrays reprocessed at 1.0 Hz | 33,439,657 | `ecf1903e1916a741fd8dd606deb80f968a5b1d76ba1ab5abce39739d009fd828` |
| ERP arrays reprocessed at 0.5 Hz | 33,459,625 | `ace0023411e0365e04953595659a0f40111d9524ffce2df9bc592a73f85d54f6` |
| ERP arrays reprocessed at 0.1 Hz | 29,014,958 | `dd14acc68cf4e6857f4ab49cfd4458ddb416be7e27a890af5dab13aa38a1b7c3` |
| `recording_qc.csv` — processing status of the high-pass reanalysis | 27,243 | `ae91856468c791952ba6d26de4f677f22e45c024fe912a2c13b246d81a6776f5` |
| `stimulus_history_individual.csv` — participant-level means by stimulus transition | 24,202 | `28bcbc711b40df10cb1a3548cd851f925a23bf46cf7fa9831b14367bec60602e` |
| `speed_order_records.csv` — walking speed and condition order | 572 | `764907b505c7610a7b35d15d179b04b3ca1fbd53ac9a805f022b416ecda10a07` |
| `raw_trigger_denominators.csv` — trigger and boundary counts | 4,735 | `fc3a795e9e3e06da6ab5f0c65829fd1183a202f575a0bc5ff794d9d1ed857d2a` |
| preprocessing log of the primary analysis | 7,984 | `e562be33d07c0b2aa8658a4253476f4e0c4c79dae3a1bd556a09df30afd2e722` |
| Event lists, `Sxx_<condition>_events.csv` (88 files) | — | SHA-256 of the 88 file hashes joined by newlines in file-name order: `136122640e6d750d04107d73c4a67f55383a1a1c3e88b62866ca77f37839880a` |

Only the canonical ERP arrays are part of the public data deposit.

## Comparison files (archived outputs)

| File | SHA-256 |
|---|---|
| `statistics/hpf_summary_python.csv` | `1ec4e10450d6274bae8c6f7bf45feb1e286c781049b2ca363b342cd763667250` |
| `statistics_converged_common/hpf_summary_python.csv` | `3c3f0de2c74a836406adeeda344cdadc60a5357a93f9a7a8b391461f53be9868` |
| `statistics/hpf_paired_filter_changes_matlab.csv` | `58d392a80cd1d5961fc4ecabe04873528ff3ca091bc99f5930a00852c64f19db` |
| `statistics_converged_common/hpf_paired_filter_changes_matlab.csv` | `1dcedf24f4f77f8a4def17111919c393987876e1b1c060ffa3000dc1654a8be1` |
| `multiplicity_audit.csv` | `2a008f1594648806e0afa3b9903b7a4eee4654484fbc78951c657387917473cd` |
| `standard_n1_report.json` | `0c7f3d83a8e828b6ffcff0e8e96bb8efae296e2a39c9b4d56169e521a5f56ccd` |
| `interaction_contrasts_report.json` | `f2188809c1a3469425f1b6dd0eee76f0d3f49ad4b3e43a99667b53704240bb22` |
| `interaction_contrasts.csv` | `01b938c3e2ab9f1dec2ae7c622dfed3f064cd9f01badc9b33a3e29d9ab145cfc` |
| `covariate_slopes_python.csv` | `05c6a49d8572a0f52e524dbdf32948379008fcf71bbae702a5dc9b25bb8a7ab9` |
| `time_window_sensitivity.csv` | `86838ddc485e0c87dc7babbc5bfdfa946f05e0c6f379fa5f5da668be4fef08ca` |
| `time_window_individual.csv` | `327a6e1e524eeac12a7acaa0a365ed52000459634b9827edc546b74e1d50175a` |
| `reviewer2_positivity_summary.csv` | `e5f1f394b7b1192bcd03db2be9c1311ec30cfdc6ffe888a50183a3b60f0bbb7d` |
| `reviewer2_positivity_individual.csv` | `360c171daf02fb147ea4017404a84749a5e67d3f1ea4b276f4f1833958e81660` |
| `early_response_descriptive.csv` | `a0b397444d2d8688f5a5f8b21879505d42d326e962d57e5e8efe7b0f4634999d` |
| `stimulus_history_report.json` | `88ff40019d66c9eb9a688c8de1f6f968111a8f2b2e31cbf608f7c5b391d058f7` |
| `speed_order_report.json` | `a056698d0373cb76b3e610d8c2c736a73ce7d5d3853095efe1db1cf0985e1c3a` |
| `sequence_rules_records.csv` | `bb85990a8e4421811189710db5fdb48fcde3c31c6eed1a0cc97489830e3e15c7` |
| `sequence_rules_report.json` | `3eb5900c4f0674e73fd40bacf385a675f8205baaab69f103e4d794b2e6c28b51` |
| `canonical_processing_participants.csv` | `a5697f93bdda17bf0710f1010331c69a142f5786918ad3fc1c1302f168465f4c` |
| `canonical_processing_summary.csv` | `fa922d5ad3ad4254a6040c0ebafe582ed450438d36330756036b73a550b704b0` |

For the reported-precision check, the comparison values were the cells of Table 3 Panel B
and Tables S2–S9 and the values given in Results 3.1, 3.2 and 3.4, the opening of Results,
Supplementary Methods S-M5 and S-M9, Supplementary Results S-R1 and the response letter.

## Tolerances and largest differences

- Estimates, confidence intervals, effect sizes, test statistics, p values and participant-level
  values: absolute tolerance 1 × 10⁻⁹; largest observed difference 3.5 × 10⁻¹⁴.
- Shapiro–Wilk p values: absolute tolerance 1 × 10⁻³ (MATLAB has no built-in function, so
  Royston's 1992 approximation was implemented); largest observed difference 3.3 × 10⁻⁹.
- Counts (stimulus-sequence totals): exact.

| Comparison | Values | Largest absolute difference | Tolerance | Archived file |
|---|---|---|---|---|
| Filter comparison: estimates, CI, dz, F, p (Table 3B, S2) | 72 | 1.6e-14 | 1e-09 | `hpf_summary_python.csv` |
| Filter comparison: Shapiro-Wilk p | 9 | 2.0e-10 | 0.001 | `hpf_summary_python.csv` |
| Filter comparison: paired changes (Table S3) | 24 | 4.4e-16 | 1e-09 | `hpf_paired_filter_changes_matlab.csv` |
| N1/MMN tests and Holm (Table S4) | 40 | 4.1e-15 | 1e-09 | `multiplicity_audit.csv` |
| N1/MMN tests: selected p and its Holm values | 12 | 4.1e-15 | 1e-09 | `multiplicity_audit.csv` |
| N1/MMN tests: Shapiro-Wilk p | 4 | 1.3e-09 | 0.001 | `multiplicity_audit.csv` |
| Standard-evoked change (Results 3.1) | 7 | 1.3e-15 | 1e-09 | `standard_n1_report.json` |
| Mean N1 interaction contrast and CI (Fig. 2e) | 4 | 4.4e-15 | 1e-09 | `interaction_contrasts_report.json` |
| Participant interaction contrasts (22 values) | 22 | 7.8e-16 | 1e-09 | `interaction_contrasts.csv` |
| Covariate slopes (Table S6) | 24 | 1.6e-14 | 1e-09 | `covariate_slopes_python.csv` |
| Latency windows (Table S7) | 140 | 5.8e-15 | 1e-09 | `time_window_sensitivity.csv` |
| Latency windows: Shapiro-Wilk p | 14 | 3.3e-09 | 0.001 | `time_window_sensitivity.csv` |
| Latency windows: participant contrasts (308 values) | 308 | 3.1e-15 | 1e-09 | `time_window_individual.csv` |
| Early positivity 40-70 ms (Table S8, S-R1) | 40 | 1.6e-15 | 1e-09 | `reviewer2_positivity_summary.csv` |
| Early positivity: participant means (88 values) | 88 | 2.2e-16 | 1e-09 | `reviewer2_positivity_individual.csv` |
| Early waveform in the common 19 (values underlying Fig. S8) | 240 | 1.5e-14 | 1e-09 | `early_response_descriptive.csv` |
| Stimulus-history summaries (Table S9) | 176 | 2.9e-14 | 1e-09 | `stimulus_history_report.json` |
| Stimulus-history paired contrasts (S-R1; response to R2 Comment 1) | 64 | 3.5e-14 | 1e-09 | `stimulus_history_report.json` |
| Walking speed and condition-by-position counts | 21 | 0.0e+00 | 1e-09 | `speed_order_report.json` |
| Stimulus sequences: per-recording counts (88 recordings) | 792 | 0.0e+00 | 1e-09 | `sequence_rules_records.csv` |
| Stimulus sequences: totals | 6 | 0.0e+00 | exact | `sequence_rules_report.json` |
| Processing counts: per recording (176 rows) | 1408 | 0.0e+00 | 1e-09 | `canonical_processing_participants.csv` |
| Processing counts: condition summaries (Table S5) | 72 | 8.9e-16 | 1e-09 | `canonical_processing_summary.csv` |

Participant-level comparisons: 22 interaction contrasts, 308 latency-window contrasts,
88 early-positivity means, sequence counts for 88 recordings, processing counts for 88 recordings.

## Reported-precision check by location

| Location | Cells | Agree |
|---|---|---|
| Table 3 Panel B | 48 | 48 |
| Table S2 | 18 | 18 |
| Table S3 | 24 | 24 |
| Table S4 | 24 | 24 |
| Table S5 and its note | 24 | 24 |
| Table S6 | 9 | 9 |
| Table S7 | 84 | 84 |
| Table S8 | 20 | 20 |
| Table S9 | 80 | 80 |
| Results 3.1 | 13 | 13 |
| Results 3.4 | 8 | 8 |
| Opening of Results (walking speed) | 4 | 4 |
| Supplementary Methods S-M5 | 8 | 8 |
| Supplementary Methods S-M9 | 4 | 4 |
| Supplementary Results S-R1 | 1 | 1 |
| Response letter | 5 | 5 |
| Total | 374 | 374 |

## Not verified here

- Table 3 Panel D (alternative pipeline without ICA). It requires ERP arrays generated from
  the raw EEG without ICA, which were not available on the machine used for this check.
- Steps that start from raw recordings: EEG preprocessing with FastICA (primary analysis and
  high-pass reanalysis), the stimulus-history preprocessing, event extraction, trigger
  counting and kinematic extraction. Their saved outputs were used as inputs.
- The transcription of walking speed and condition order from the laboratory record sheet
  into `speed_order_records.csv`.
- Drawing of the figures.
