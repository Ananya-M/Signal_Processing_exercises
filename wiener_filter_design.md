# Wiener filter design for SSEP analysis (SOMA dataset)

## Overview

This document captures the design rationale, implementation decisions, pitfalls encountered, and fixes applied when building a frequency-domain Wiener filter for somatosensory evoked potential (SSEP) analysis. It is intended as a reference for anyone diving into adaptive filtering replicating the pipeline.

---

## 1. Theory recap

### The Wiener filter in one line

The optimal linear filter that minimises mean squared error between a desired signal $d[n]$ and a noisy observation $x[n]$ is:

$$H^*(f) = \frac{S_{dd}(f)}{S_{dd}(f) + S_{vv}(f)}$$

where $S_{dd}$ is the power spectral density of the signal and $S_{vv}$ is the PSD of the noise. This is the **Wiener-Khinchin** form: the filter gain at each frequency equals the signal-to-noise ratio normalised to [0, 1].

- Where signal dominates: $H^* \to 1$ — pass through unchanged
- Where noise dominates: $H^* \to 0$ — suppress
- The filter is data-driven: its shape is determined entirely by the spectral content of your data

### Why this suits SSEP

SSEPs are phase-locked to the stimulus and repeat across trials. The cross-trial average is a clean estimate of the deterministic SSEP template. Trial-to-trial residuals are dominated by ongoing EEG noise (alpha, beta, drift). The Wiener filter exploits exactly this structure.

### The Wiener-Hopf equation (time domain equivalent)

$$\mathbf{R}_{xx}\,\mathbf{h}^* = \mathbf{r}_{dx}$$

Working in the frequency domain avoids solving this linear system explicitly and is $O(T \log T)$ via the FFT rather than $O(M^3)$ for matrix inversion.

---

## 2. Pipeline design

### Inputs

| Variable | Shape | Description |
|---|---|---|
| `LB` | `[T × N]` | Local baseline trials — T time samples (rows), N trials (columns) |
| `Fs` | scalar | Sampling rate (Hz) |

**Critical dimension rule:** rows = time, columns = trials. If T < N you have a transposed matrix. The function checks this and warns.

### Steps

```
1. Template d[n]    = nanmean(LB, 2)          [T×1]   SSEP estimate
2. Residuals V      = LB − d                  [T×N]   noise estimate
3. Shared Hann win  = hann(T)                 [T×1]   same for both PSDs
4. Sdd              = |FFT(d .* win)|² / scale        SSEP PSD
5. Sdd              = Sdd * N                          scale to per-trial power
6. Svv              = mean over trials of |FFT(v .* win)|² / scale
7. H_raw            = Sdd ./ (Sdd + Svv + eps)        Wiener gain
8. H_smooth         = Gaussian smooth of H_raw         stabilise ratio
9. taper            = cosine rolloff at band edges     remove ringing
10. H               = H_smooth .* taper                final gain
11. h               = real(ifft(H))                    impulse response
```

### Applying the filter

```matlab
X_fft  = fft(X_clean, [], 1);      % FFT along time dimension
Y_fft  = X_fft .* H;               % element-wise multiply (broadcast)
Y_time = real(ifft(Y_fft, [], 1)); % back to time domain
```

Multiplying in the frequency domain is equivalent to circular convolution with `h = ifft(H)` (Convolution Theorem). For short epochs with a smooth H, circular convolution is a safe approximation.

---

## 3. Pitfalls encountered and fixes applied

### Pitfall 1 — Dimension swap (T and N transposed)

**What happened:** The function assigned `N = size(LB, 1)` and `T = size(LB, 2)`, swapping time samples and trial count. The Hann window had length N (number of trials) instead of T (number of samples), and `nanmean(LB, 2)` averaged across time instead of across trials.

**Symptom:** Everything downstream was wrong — PSDs were garbage, the filter had no meaningful shape, filtered ERPs were flat or inverted.

**Fix:**
```matlab
T = size(LB, 1);   % time samples — rows
N = size(LB, 2);   % trials       — columns
if T < N
    warning('Matrix may be transposed: T=%d, N=%d', T, N);
end
```

**Rule:** Always label dimension assignments with a comment. T should be much larger than N for EEG epochs.

---

### Pitfall 2 — Epoch too long, signal diluted

**What happened:** The filter was estimated over the full −400 to +800 ms epoch (~1200 ms). The SSEP (N20/P35 complex) occupies only the first ~80 ms. For the remaining ~1100 ms, `d[n]` is near zero. The periodogram of `d` was therefore tiny relative to `S_vv` (broadband EEG noise present throughout the full epoch), making `H* ≈ 0` everywhere.

**Symptom:** H*(f) never exceeded 0.1 anywhere. Filtered ERPs were flat at zero.

**Fix:** Trim to the SSEP-relevant window before calling the filter:
```matlab
t_win = EEG.times >= 10 & EEG.times <= 125;   % ms
LB_ssep = LB(t_win, :);
PS_ssep = PS(t_win, :);
[H_w, h_w, Sdd, Svv, f_axis] = compute_wiener_filter(LB_ssep, Fs);
```

**Rule:** The epoch passed to the Wiener filter should be restricted to the time window where signal actually exists. A longer epoch is not better — it dilutes the signal PSD.

---

### Pitfall 3 — Sdd not scaled to per-trial power

**What happened:** `S_dd` is the PSD of the cross-trial average `d = mean(LB, 2)`. Averaging N trials reduces the amplitude of the incoherent noise by $1/\sqrt{N}$ and the power by $1/N$. So `S_dd` is N times smaller than the per-trial signal power. Meanwhile `S_vv` is estimated per trial and averaged — it is on the per-trial scale. The ratio `S_dd / (S_dd + S_vv)` was therefore always near zero.

**Symptom:** H*(f) near 0 everywhere even after epoch trimming. Both PSDs appeared tiny in the diagnostic plot.

**Fix:**
```matlab
Sdd = (abs(fft(d .* win)).^2) / scale;
Sdd = Sdd * N;   % restore per-trial scale
```

**Rule:** `S_dd` must be scaled up by N to be comparable to the per-trial `S_vv`. Without this the SNR estimate is always artificially low by a factor of N.

---

### Pitfall 4 — Asymmetric smoothing of Sdd and Svv

**What happened:** `S_dd` was smoothed with `movmean(Sdd, 15)` (wide) and `S_vv` with `movmean(Svv, 5)` (narrow). The ratio then inherited the shape of $1/S_{vv}$ — a high-pass that rises wherever noise falls — rather than reflecting true SNR. The filter looked physiological but was driven by smoothing mismatch, not signal content.

**Symptom:** H*(f) monotonically increased from 0 to 1 — a high-pass shape rather than the expected band-pass dome.

**Fix:** Smooth the ratio H directly, not the inputs separately:
```matlab
H_raw = Sdd ./ (Sdd + Svv + eps);

sigma = 3;
x_k   = (-3*sigma : 3*sigma)';
gauss = exp(-x_k.^2 / (2*sigma^2));
gauss = gauss / sum(gauss);
H_smooth = conv(H_raw, gauss, 'same');
H_smooth = max(H_smooth, 0);
```

**Rule:** Never smooth inputs asymmetrically. If smoothing is needed, apply it to the ratio. Gaussian kernel preferred over `movmean` — softer rolloff, no boxcar edges.

---

### Pitfall 5 — Normalisation mismatch between Sdd and Svv

**What happened:** `S_dd` used one normalisation (no Hann window, divide by T) and `S_vv` used another (Hann window, divide by `T * win_norm`). Even though the scale cancels in the ratio, the mismatch changed the effective relative magnitudes, artificially suppressing `S_dd` relative to `S_vv`.

**Symptom:** H* driven to 1 at high frequencies where both PSDs were near zero, because `eps` dominated the denominator — a numerical artefact, not a real filter.

**Fix:** Use the same window and the same scale for both:
```matlab
win   = hann(T);
scale = sum(win.^2);   % identical for both

D   = fft(d .* win);
Sdd = (abs(D).^2) / scale;   % consistent

V_n = fft(v_n .* win);
Svv = Svv + (abs(V_n).^2) / scale;   % same scale
```

**Rule:** The scale cancels in the ratio, so absolute units don't matter. But relative normalisation does. Use the same window and divisor for `S_dd` and `S_vv`.

---

### Pitfall 6 — Gibbs ringing from the preprocessing brick wall

**What happened:** The preprocessing applied a 50 Hz lowpass, so H(f) drops abruptly to zero at 50 Hz. The `ifft` of a sharp step in frequency is a sinc — infinite ringing in time. This appeared as a 20 ms period sinusoid spanning the full ±60 ms impulse response window (period = 1/50 Hz = 20 ms).

**Symptom:** Panel 4 (impulse response) showed a wide sinusoidal oscillation. Filtered ERPs had artefactual ringing that obscured the N20/P35 complex.

**Fix:** Apply a cosine taper to smooth H to zero before the cutoff:
```matlab
f_high_start = 35;   % Hz — begin rolloff
f_high_zero  = 50;   % Hz — reach zero (matches preprocessing cutoff)

idx_hi = f_axis >= f_high_start & f_axis <= f_high_zero;
taper(idx_hi) = 0.5 * (1 + cos(pi * ...
    (f_axis(idx_hi) - f_high_start) / (f_high_zero - f_high_start)));
```

The formula gives taper = 1 at 35 Hz, 0.5 at 42 Hz, and 0 at 50 Hz — a smooth S-curve with no discontinuity.

**Rule:** Never let H drop abruptly to zero at the band edge. Always apply a smooth taper (cosine, Tukey, or Gaussian) over a transition zone of at least 10–15 Hz.

---

### Pitfall 7 — Taper not mirrored to negative frequencies

**What happened:** The cosine taper was applied only to the positive frequency side (35–50 Hz). The FFT spectrum repeats after Fs/2 — the bins from Fs−50 to Fs−35 Hz are the negative-frequency mirror of the 35–50 Hz rolloff. Without mirroring, H had a smooth rolloff on the positive side but a hard jump on the negative side, breaking conjugate symmetry.

**Consequence:** `ifft(H)` had a non-zero imaginary part. Taking `real(ifft(H))` discarded the imaginary part incorrectly, producing an asymmetric impulse response (left lobe ≠ right lobe).

**Fix:** Mirror the taper to both sides of the spectrum:
```matlab
% positive side: 35 → 50 Hz (taper down)
idx_hi = f_axis >= f_high_start & f_axis <= f_high_zero;
taper(idx_hi) = 0.5 * (1 + cos(pi * ...
    (f_axis(idx_hi) - f_high_start) / (f_high_zero - f_high_start)));

% hard zero in the middle
idx_zero = f_axis > f_high_zero & f_axis < (Fs - f_high_zero);
taper(idx_zero) = 0;

% negative side: Fs−50 → Fs−35 Hz (taper up — mirror)
idx_neg = f_axis >= (Fs - f_high_zero) & f_axis <= (Fs - f_high_start);
taper(idx_neg) = 0.5 * (1 + cos(pi * ...
    ((Fs - f_axis(idx_neg)) - f_high_start) / (f_high_zero - f_high_start)));
```

**Rule:** Any modification to H on the positive frequency side must be mirrored symmetrically on the negative frequency side. The FFT of a real signal satisfies `H[k] = conj(H[T-k])`. Violating this makes `ifft(H)` complex.

---

### Pitfall 8 — Conflicting hard cutoffs overwriting the taper

**What happened:** The code applied a cosine taper from 35–48 Hz and then separately zeroed all bins above 45 Hz. The hard zero at 45 Hz overwrote the taper for bins between 45 and 48 Hz — the taper in that region was computed and then immediately discarded.

**Symptom:** Effective cutoff was at 45 Hz with a partial taper from 35–45 Hz, plus an unintended hard edge at 45 Hz.

**Fix:** One taper, one cutoff. The taper ends exactly at the cutoff; the hard zero starts exactly where the taper ends:

```matlab
f_high_start = 35;   % taper begins
f_high_zero  = 50;   % taper reaches zero AND hard zero begins

% These two regions are contiguous and non-overlapping:
% [35, 50]   → cosine taper
% (50, Fs-50) → hard zero
```

**Rule:** Define `f_high_start` and `f_high_zero` once. The hard zero region must begin exactly where the taper ends — no gap, no overlap.

---

## 4. Constraints imposed by preprocessing

The SOMA dataset was preprocessed with a **1–50 Hz bandpass** before epoching. This has two important consequences:

1. **The Wiener filter cannot recover frequencies above 50 Hz.** The N20/P35 complex has energy up to ~150 Hz in wideband recordings. After a 50 Hz lowpass, only the slow envelope of the SSEP survives. The filter operates on this broadened shape, not the sharp N20 peak.

2. **The useful frequency resolution is limited.** With a short epoch (~115 ms) and a 50 Hz bandwidth, you have roughly `50 / (Fs/T)` ≈ 6–12 meaningful frequency bins. Smoothing the ratio H is essential to avoid noise in the gain estimate from so few bins.

**Implication for metrics:** Peak amplitude and latency of the N20 are less precise than in wideband recordings. Area under the curve of the N20/P35 window (15–45 ms) is a more robust metric in this bandwidth.

**Implication for the filter shape:** A correctly working filter on this data should show H*(f) with a dome shape peaking in the SSEP band. If H*(f) is monotonically increasing (high-pass shape), the filter is being driven by noise statistics rather than signal content — a sign that either the epoch is too long or `S_dd` is not scaled correctly.

---

## 5. Diagnostic checklist

After each run, inspect the four panels:

| Panel | What to look for | Red flag |
|---|---|---|
| Panel 1 (ERPs) | Wiener ERP tracks raw ERP in SSEP band, suppresses noisy tail | Wiener ERP flat, inverted, or oscillating |
| Panel 2 (H*(f)) | Dome shape peaking in SSEP band, smooth rolloff at 35–50 Hz | Monotonic high-pass, flat at 1, or erratic |
| Panel 3 (PSDs) | S_dd visibly above S_vv in SSEP band | Both curves invisible or S_vv always above S_dd |
| Panel 4 (h[n]) | Compact single lobe centred at lag 0, symmetric, decays by ±20 ms | Sinusoidal ringing, asymmetric lobes, wide spread |

**Residual power:** BS and PS residual powers should differ if LIFU has an effect. If they are equal, the filter is not recovering signal — check upstream issues.

---

## 6. Key parameters

| Parameter | Value | Rationale |
|---|---|---|
| Epoch window | 10–125 ms | SSEP-relevant window only; avoids diluting S_dd |
| Window function | Hann | Reduces spectral leakage; same for S_dd and S_vv |
| Sdd scaling | × N trials | Restores per-trial power after cross-trial averaging |
| Smoothing kernel | Gaussian, σ = 3 bins | Applied to ratio H directly; softer than movmean |
| Taper start | 35 Hz | 15 Hz transition zone before preprocessing cutoff |
| Taper end / hard zero | 50 Hz | Matches preprocessing bandpass |
| DC cutoff | 2 Hz | Removes residual drift |
| Mirror taper | Fs−50 to Fs−35 Hz | Restores conjugate symmetry of H |

---

## 7. What this filter cannot do

- **Recover frequencies removed by preprocessing.** If the data was lowpassed at 50 Hz, no filter can restore the 50–150 Hz content of the N20.
- **Handle non-stationary noise.** The filter is estimated from baseline and applied to post-sonication. If noise statistics change (e.g. increased muscle artefact post-LIFU), the filter is no longer optimal.
- **Distinguish LIFU effects from noise.** If post-sonication SNR drops to near zero, the filter output is shaped noise, not a filtered SSEP. Always check residual power and compare to baseline.
- **Replace trial rejection.** Heavily artefacted trials should be removed before the filter is estimated. Including them raises S_vv artificially and suppresses the filter gain.
