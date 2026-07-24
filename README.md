# walking-eeg-n1-mmn

MATLAB re-implementation of the statistical analyses reported in the manuscript *"Preferential attenuation of deviant-evoked auditory N1 responses during treadmill walking"* (auditory N1/MMN, passive frequency oddball; EEG n = 22, kinematics n = 16–17). The scripts reproduce specified tables, panels, or rows deterministically from the deposited reproducibility data. The Table 4 Panel B and Panel D sensitivity analyses additionally require ERP-array variants available from the corresponding author; see "Scope and two on-request items" below.

This repository contains **code only**. The data are deposited separately and will be made openly available upon publication (during peer review they are available to editors and reviewers via a private link); see the Data section below.

## Data

The data underlying the tables reported in the manuscript will be made openly available on OpenNeuro **ds008198**
(https://openneuro.org/datasets/ds008198) upon publication; during peer review they are available
to editors and reviewers via a private link. The deposit comprises raw EEG (BIDS-EEG), the raw wrist/L3
accelerometry (`sourcedata/kinematic/`), and a reproducibility set
(`derivatives/reproducibility/`) with the subject-/condition-averaged ERP arrays,
the derived kinematic outcome tables, and the retained trial counts. The Table 4 Panel B and
Panel D sensitivity analyses additionally require ERP-array variants available from the
corresponding author (see "Scope and two on-request items" below).

Put these deposited files (from `derivatives/reproducibility/` on OpenNeuro) in a
`data/` folder next to `code/` (the scripts default to `../data`):

| File | Backs |
|---|---|
| `erp_arrays_n22.mat` | subject × condition × channel × time ERP arrays (separate standard-, deviant-, and difference-evoked arrays) — all reported EEG results except the Table 4 Panel B and Panel D variant rows |
| `kinematic_wrist_3cond_paired_wide_n22.csv` | Table 3B (wrist) |
| `kinematic_L3_3cond_paired_wide_n22.csv` | Table 3B (L3) |
| `per_subject_trial_counts_n22.csv` | Table 1B, Table 4 Panel C |

**Determinism.** Every statistic reproduced here is deterministic *given the deposited
derived data* (the two Table 4 sensitivity variants, Panels B and D, instead require the
on-request ERP-array variants; see below). Two upstream steps are tool/platform-dependent and are therefore
deposited as derived tables rather than regenerated here: the raw→ERP-arrays step
(FastICA) and the raw→derived-kinematics step (drift-free Fourier displacement with
peak-based step detection). The scope of this repository is statistical reproduction from the
deposited derived data; regeneration of the ERP arrays and kinematic outcome tables from the
raw recordings is not included. The canonical anchor is the N1 Sit × Walk-Free interaction
**F(1,21) = 6.26** (partial η² = .230).

## Requirements

MATLAB R2025b (Statistics and Machine Learning Toolbox: `fitlme`, `signrank`, `tcdf`, …).

## Scripts

Run each as `fname('<path>')` (or with no argument if the data are in `../data`):

```
matlab -batch "stats_analysis1_sit_walk('data/erp_arrays_n22.mat')"
```

| Script | Reproduces | Argument(s) |
|---|---|---|
| `stats_analysis1_sit_walk.m` | **Table 2** — Sit vs Walk-Free N1-window analysis | ERP `.mat` |
| `stats_analysis1_n1_interaction.m` | Table 2 — N1 Stimulus × Condition interaction | ERP `.mat` |
| `stats_analysis2_water_clay.m` | **Table 3A** (MMN amplitude + TOST) and **Supp. Table S1** (permutation/Wilcoxon corroboration) | ERP `.mat` |
| `check_topo_consistency.m` | **Supp. Table S1** (topographic-consistency rows) | ERP `.mat` |
| `stats_table3a_gmd_tanova.m` | **Supp. Table S1** (GMD / exact-TANOVA scalp-map row) | ERP `.mat` |
| `stats_kinematic_table3b.m` | **Table 3B** — Water vs Clay kinematics (4 outcomes) | `data/` dir |
| `stats_table1b_trialcounts.m` | **Table 1B** — retained trial-count RM-ANOVA | `data/` dir |
| `stats_table4a_timewindow.m` | **Table 4 Panel A** — time-window sensitivity | ERP `.mat` |
| `stats_table4c_covariate.m` | **Table 4 Panel C** — trial-count covariate check | ERP `.mat`, counts `.csv` |

Each script prints its results next to the target values reported in the manuscript in `[...]`.

## Methods notes (why the numbers match to the last digit)

- **Wilcoxon signed-rank** uses the *exact* distribution: `signrank(d, 0, 'method', 'exact')`.
  MATLAB's default `signrank` switches to a normal approximation for n > 15, which does
  not reproduce the reported (exact) p-values.
- **Permutation p** (the sign-flip corroboration tests and the label-swap TANOVA) is the
  *exact* test — all 2^n reassignments are enumerated, so the value is deterministic and
  platform-independent (no RNG seed).
- **Kinematics.** The L3 excursion inputs are orientation-corrected, drift-free per-stride Fourier
  displacements. No additional exclusions based on the derived kinematic values are applied by the
  scripts; the paired samples are wrist n = 16 and L3 n = 17, as determined by sensor-data availability.

## Scope and two on-request items

- `stats_table4c_covariate.m` reproduces Table 4 Panel C using a transparent
  within-participant contrast ANCOVA (base β = +0.280, F(1,21) = 6.26; covariate-adjusted
  β = +0.280, F(1,20) = 6.78). The interaction remained significant after adjustment for
  retained trial counts.
- **Table 4 Panel B** (0.5- and 0.1-Hz high-pass filters) and **Panel D** (no ICA) require the
  corresponding ERP-array variants generated by the original preprocessing pipeline. These
  variants are available from the corresponding author upon reasonable request. The
  1.0-Hz/FastICA reference rows can be reproduced from the deposited canonical ERP array.
