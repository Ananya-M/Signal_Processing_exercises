# Signal Processing Exercises

Hands-on MATLAB exercises in signal processing for EEG, built around somatosensory evoked potentials (SSEPs). Each exercise has its own folder with a script and a README covering the theory, how the pipeline is built, common mistakes and how to read the results.

## Exercises

| # | Exercise | Topics | Status |
|---|---|---|---|
| 01 | [Wiener filter](exercises/01_wiener_filter/) | Frequency-domain Wiener filter, PSD estimation, Wiener-Khinchin, spectral tapering, impulse response diagnostics | Available |
| 02 | Adaptive filtering (LMS) | Removing ECG artifacts from continuous EEG with NLMS | Coming soon |
| 03 | Surface Laplacian | Spherical-spline surface Laplacian vs. conventional re-referencing | Coming soon |

## Repository structure

```
Signal_Processing_exercises/
├── README.md                  ← you are here
└── exercises/
    └── 01_wiener_filter/
        ├── README.md          ← theory, pitfalls, diagnostic checklist
        ├── wiener_filter_ssep.m
        └── images/
```

## Getting started

1. Install MATLAB (with the Signal Processing Toolbox) and [EEGLAB](https://sccn.ucsd.edu/eeglab/).
2. Clone this repository:
   ```
   git clone https://github.com/Ananya-M/Signal_Processing_exercises.git
   ```
3. Choose an exercise and read its README before running the script.

## Data

The exercises use a sample EEG dataset that isn't stored in this repository. If you'd like to try the exercises, email me and I'll send it to you. Each exercise README explains where to point the script once you have the data.
