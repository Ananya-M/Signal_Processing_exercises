# EEG Signal Processing to try out things that I didn't in my pipeline?

MATLAB exercises(for now) in signal processing for EEG, built around somatosensory evoked potentials (SSEPs). Each exercise has its own folder with a script and a README covering the theory, how the script is built, some outputs that I got and what I think about them.

## Exercises

| # | Exercise | Topics | Status |
|---|---|---|---|
| 01 | [Wiener filter](exercises/01_wiener_filter/) | Frequency-domain Wiener filter, PSD estimation, Wiener-Khinchin, spectral tapering, impulse response diagnostics | Available |
| 02 | [Adaptive filtering (LMS)](exercises/02_adaptive_filtering_lms/) | Removing ECG artifacts from continuous EEG with NLMS, R-peak/stimulus phase-locking (Rayleigh test) | Available |
| 03 | Surface Laplacian | Spherical-spline surface Laplacian vs. conventional re-referencing | Coming soon |

## Repository structure

```
Signal_Processing_exercises/
├── README.md                  ← you are here
├── setup_paths.m              ← adds utils/ to the MATLAB path
├── utils/                     ← shared helper functions
└── exercises/
    ├── 01_wiener_filter/
    │   ├── README.md          ← theory, pitfalls, diagnostic checklist
    │   ├── wiener_filter_ssep.m
    │   └── images/
    └── 02_adaptive_filtering_lms/
        ├── README.md          ← LMS/NLMS theory, filter design, phase-locking check
        ├── nlms_ekg_removal.m
        └── images/
```

## Getting started

1. Install MATLAB (with the Signal Processing Toolbox) and [EEGLAB](https://sccn.ucsd.edu/eeglab/).
2. Clone this repository:
   ```
   git clone https://github.com/Ananya-M/Signal_Processing_exercises.git
   ```
3. In MATLAB, run `setup_paths` from the repository root.
4. Choose an exercise and read its README before running the script.

## Data

The exercises use a sample EEG dataset that isn't stored in this repository. If you'd like to try the exercises, email me and I'll send it to you. Each exercise README explains where to point the script once you have the data.
