clc; clear;
%% ========== MAIN LOOP ==================================================


tic;
subject_id="Sample_data";
SUBJECT = struct();

data_dir = 'path/to/Sample_data';   % TODO: folder containing Sample_data.set and Sample_data.fdt
EEG = pop_loadset('Sample_data.set', data_dir);
EEG = eeg_checkset(EEG);
Fs  = EEG.srate;
first_protocol=fieldnames(EEG.protocoldata.Trigger_13);
first_protocol=first_protocol{1};
fp=matlab.lang.makeValidName(first_protocol);

time=EEG.times;
idx_10  = round(interp1(time, 1:length(time), 10));
idx_175 = round(interp1(time, 1:length(time), 175));

EEG.protocoldata.Trigger_13.Global_pre.baseline=EEG.protocoldata.Trigger_13.(fp).pre;
EEG.protocoldata.Trigger_14.Global_pre.baseline=EEG.protocoldata.Trigger_14.(fp).pre;
protocols = fieldnames(EEG.protocoldata.Trigger_13);
protocols = ["Global_pre"; protocols(1:end)];
protocols(end) = [];

sides = {'Trigger_13', 'Trigger_14'};

for si = 1%:numel(sides)
    side = sides{si};

    if strcmp(side, 'Trigger_13')
        group    = 'Task';
        channels = {'C4','CP2','CP6'};
    else
        group    = 'CONTROL';
        channels = {'C3','CP1','CP5'};
    end

    for p = 2%:numel(protocols)
        prot = protocols{p};

        for ci = 3 %1:numel(channels)
            chan = channels{ci};

            %% --- 1. RAW EPOCH MATRICES ----------------------------
            GB = EEG.protocoldata.(side).(prot).pre.(chan);
            if size(GB, 2) > 300
                GB = GB(idx_10:idx_175, end-200:end);
            end

            LB = EEG.protocoldata.(side).(prot).pre.(chan);

            if p==2
                LB = LB(idx_10:idx_175, end-110:end);
            else
                LB = LB(idx_10:idx_175, :);
            end

            % Before calling compute_wiener_filter, add:
            figure;
            subplot(2,1,1); imagesc(LB'); colorbar; title('local pre raw trials [N x T]');
            subplot(2,1,2); plot(nanmean(LB,2)); title('ERP template d[n]');

            %% --- 2. CONVENTIONAL ERPs -----------------------------
            ERP_gb_uV = nanmean(GB, 2);
            ERP_lb_uV = nanmean(LB, 2);


            %% --- 3. WIENER FILTER (fully frequency-domain) --------
            [H_w, h_w, Sdd, Svv, f_axis] = compute_wiener_filter(LB, Fs,time(idx_10:idx_175));

            LB_wiener = apply_wiener_epoch(LB, H_w);
            ERP_lb_w = nanmean(LB_wiener, 2);
            resid_power_lb = mean((LB_wiener - ERP_lb_w).^2, 'all', 'omitnan');


            if ~strcmp(prot,'Global_pre')

                PS = EEG.protocoldata.(side).(prot).post.(chan);
                PS = PS(idx_10:idx_175, :);
                ERP_ps_uV = nanmean(PS, 2);
                PS_wiener = apply_wiener_epoch(PS, H_w);
                ERP_ps_w = nanmean(PS_wiener, 2);
                resid_power_ps = mean((PS_wiener - ERP_ps_w).^2, 'all', 'omitnan');
            end

            %% --- 4. STORE -----------------------------------------
            tag = sprintf('%s__%s__%s__%s', subject_id, group, prot, chan);

            SUBJECT.(tag).subject_id     = subject_id;
            SUBJECT.(tag).group          = group;
            SUBJECT.(tag).protocol       = prot;
            SUBJECT.(tag).channel        = chan;
            SUBJECT.(tag).ERP_gb_uV      = ERP_gb_uV;
            SUBJECT.(tag).ERP_lb_uV      = ERP_lb_uV;
            SUBJECT.(tag).ERP_lb_wiener  = ERP_lb_w;
            SUBJECT.(tag).wiener_H       = H_w;
            SUBJECT.(tag).wiener_h       = h_w;
            SUBJECT.(tag).wiener_Sdd     = Sdd;
            SUBJECT.(tag).wiener_Svv     = Svv;
            SUBJECT.(tag).wiener_f       = f_axis;
            SUBJECT.(tag).resid_power_lb = resid_power_lb;


            if ~strcmp(prot,'Global_pre')
                SUBJECT.(tag).ERP_ps_uV      = ERP_ps_uV;
                SUBJECT.(tag).ERP_ps_wiener  = ERP_ps_w;
                SUBJECT.(tag).resid_power_ps = resid_power_ps;

                if strcmp(group, "Task") && strcmp(prot, "ProtocolC")|| strcmp(prot, "ProtocolR") && strcmp(chan, "CP6")
                    plot_wiener_summary(EEG.times(idx_10:idx_175), ...
                        ERP_lb_uV, ERP_ps_uV, ...
                        ERP_lb_w,  ERP_ps_w, ...
                        h_w, Sdd, Svv, f_axis, Fs, ...
                        resid_power_lb, resid_power_ps, tag);
                end
            end

        end % channel
    end % protocol
end % side


toc;
fprintf('\nDone.\n');


%% ========================================================================
%%  FUNCTION: compute_wiener_filter
%% ========================================================================
function [H, h, Sdd, Svv, f_axis] = compute_wiener_filter(LB, Fs, time)
%COMPUTE_WIENER_FILTER  Frequency-domain Wiener filter from baseline trials.
%
%  LB     : [T x N]  local baseline — T time samples (rows), N trials (cols)
%  Fs     : sampling rate (Hz)
%  time   : [T x 1]  time axis (ms) — passed in but used only for reference
%
%  Returns
%  H      : [T x 1]  Wiener gain spectrum, real, tapered, in [0,1]
%  h      : [T x 1]  impulse response = real(ifft(H))  [diagnostic]
%  Sdd    : [T x 1]  SSEP template PSD  (per-trial scaled)
%  Svv    : [T x 1]  noise residual PSD (averaged across trials)
%  f_axis : [T x 1]  frequency axis (Hz), 0 .. Fs

    %% ── Dimensions ───────────────────────────────────────────────────────
    T = size(LB, 1);   % number of time samples  (rows)
    N = size(LB, 2);   % number of trials        (columns)

    if T < N
        warning('compute_wiener_filter: T=%d < N=%d — matrix may be transposed.', T, N);
    end

    %% ── Step 1: SSEP template and noise residuals ────────────────────────
    % Template = cross-trial average → cancels incoherent noise,
    % retains only the phase-locked SSEP component.
    d = nanmean(LB, 2);   % [T x 1]

    % Residuals = each trial minus the template → pure noise estimate.
    % NaN-safe: nanmean ignores NaNs, so subtraction is fine column-wise.
    V = LB - d;            % [T x N]

    %% ── Step 2: Shared window and scale ─────────────────────────────────
    % Hann window applied to BOTH d and residuals — identical normalisation
    % so the scale factor cancels exactly in the ratio H = Sdd/(Sdd+Svv).
    % We therefore do not need physical units (µV²/Hz); shape is all that
    % matters for the gain.
    win   = hann(T);          % [T x 1]
    scale = sum(win.^2);      % scalar — power normalisation (cancels in ratio)

    %% ── Step 3: SSEP PSD ────────────────────────────────────────────────
    % Single periodogram of the ERP average.
    % The average reduces amplitude by 1/N relative to one trial, so we
    % scale back up by N to put Sdd on the same per-trial power scale as Svv.
    % Without this, Sdd << Svv everywhere and H → 0 everywhere.
    D   = fft(d .* win);             % [T x 1]
    Sdd = (abs(D).^2) / scale * N;  % [T x 1]  per-trial scale

    %% ── Step 4: Noise PSD ───────────────────────────────────────────────
    % Bartlett estimate: average periodogram across all valid residual trials.
    % Averaging reduces variance of the PSD estimate (more trials = smoother).
    Svv     = zeros(T, 1);
    n_valid = 0;
    for n = 1:N
        v_n = V(:, n);
        if all(isnan(v_n)), continue; end
        v_n(isnan(v_n)) = 0;          % replace isolated NaNs with 0
        V_n  = fft(v_n .* win);       % [T x 1]
        Svv  = Svv + (abs(V_n).^2) / scale;
        n_valid = n_valid + 1;
    end
    if n_valid > 0
        Svv = Svv / n_valid;          % [T x 1]
    end

    %% ── Step 5: Frequency axis ───────────────────────────────────────────
    f_axis = (0:T-1)' * (Fs / T);    % [T x 1], runs 0 → Fs in steps of Fs/T

    %% ── Step 6: Raw Wiener gain ──────────────────────────────────────────
    % H*(f) = Sdd / (Sdd + Svv)
    % Both Sdd and Svv are real and ≥ 0.
    % eps prevents 0/0 at bins where both are zero (above preprocessing cutoff).
    H_raw = Sdd ./ (Sdd + Svv + eps);   % [T x 1], real, in [0,1]
    
    % Add this right after computing Sdd and Svv in compute_wiener_filter
    f_10hz = find(f_axis >= 10, 1);
    f_30hz = find(f_axis >= 30, 1);
    fprintf('Sdd at 10Hz: %.4f   at 30Hz: %.4f\n', Sdd(f_10hz), Sdd(f_30hz));
    fprintf('Svv at 10Hz: %.4f   at 30Hz: %.4f\n', Svv(f_10hz), Svv(f_30hz));
    fprintf('N trials: %d   T samples: %d   Fs: %.1f\n', N, T, Fs);
    %% ── Step 7: Smooth H directly ────────────────────────────────────────
    % Smoothing the ratio rather than Sdd and Svv separately avoids
    % asymmetric smoothing artefacts that create false filter shapes.
    % Gaussian kernel — softer rolloff than movmean, no ringing.
    sigma    = 3;                              % bins (~Fs/T * 3 Hz per bin)
    x_k      = (-3*sigma : 3*sigma)';
    gauss    = exp(-x_k.^2 / (2*sigma^2));
    gauss    = gauss / sum(gauss);
    H_smooth = conv(H_raw, gauss, 'same');    % [T x 1]
    H_smooth = max(H_smooth, 0);              % clip numerical negatives

    %% ── Step 8: Cosine taper ─────────────────────────────────────────────
    % Tapers H smoothly to zero at the preprocessing bandpass edges.
    % Avoids the sharp spectral discontinuity that causes Gibbs ringing
    % in the impulse response
    %
    % Low-frequency taper: DC to f_low_zero (removes drift)
    % High-frequency taper: f_high_start → f_high_zero (matches 50 Hz LP)
    % Both tapers are mirrored to the negative-frequency side.

    f_low_zero    = 2;    % Hz — hard zero below this
    f_high_start  = 35;   % Hz — taper begins rolling off here
    f_high_zero   = 50;   % Hz — reaches zero here (your preprocessing cutoff)

    taper = ones(T, 1);

    % %% Low-frequency: hard zero below f_low_zero
    taper(f_axis < f_low_zero) = 0;
    % Mirror: negative frequencies above Fs - f_low_zero
    taper(f_axis > (Fs - f_low_zero)) = 0;

    % %% High-frequency cosine taper — positive side
    idx_hi = f_axis >= f_high_start & f_axis <= f_high_zero;
    taper(idx_hi) = 0.5 * (1 + cos(pi * ...
        (f_axis(idx_hi) - f_high_start) / (f_high_zero - f_high_start)));
    % 
    % %% Hard zero between f_high_zero and Fs - f_high_zero (both sides zeroed)
    idx_zero = f_axis > f_high_zero & f_axis < (Fs - f_high_zero);
    taper(idx_zero) = 0;
    % 
    % %% High-frequency cosine taper — negative side (mirror)
    idx_neg = f_axis >= (Fs - f_high_zero) & f_axis <= (Fs - f_high_start);
    taper(idx_neg) = 0.5 * (1 + cos(pi * ...
        ((Fs - f_axis(idx_neg)) - f_high_start) / (f_high_zero - f_high_start)));

    %% Apply taper
    H = H_smooth .* taper;   % [T x 1], real, in [0,1]

    %% ── Step 9: Impulse response (diagnostic) ────────────────────────────
    % H is real and conjugate-symmetric (by construction above).
    % ifft of a real symmetric spectrum is real.
    % fftshift(h) centres the zero-phase response for plotting.
    h = real(ifft(H));        % [T x 1]

end
%% ========================================================================
%%  FUNCTION: apply_wiener_epoch
%% ========================================================================
function Y = apply_wiener_epoch(X, H)
%APPLY_WIENER_EPOCH  Multiply each trial's spectrum by H(f), return to time.
%
%  Because we work in the frequency domain the operation is:
%       Y(:,n) = ifft( fft(X(:,n)) .* H )
%
%  This is equivalent to circular convolution with h = ifft(H).
%  For SSEP epochs (short, ~100 ms) circular ≈ linear convolution
%  because the filter h is smooth and decays well within the epoch.
%
%  Inputs
%   X : [T x N]  raw EEG trials (may contain NaNs)
%   H : [T x 1]  Wiener gain spectrum (output of compute_wiener_filter)
%
%  Output
%   Y : [T x N]  filtered trials, real-valued

    [T, N] = size(X);
    Y = nan(T, N);

    % FFT of all trials at once — MATLAB fft operates column-wise by default
    % Replace NaNs with 0 before FFT (avoids NaN propagation through FFT)
    X_clean = X;
    X_clean(isnan(X)) = 0;

    X_fft = fft(X_clean, [], 1);          % [T x N] — FFT along time (dim 1)

    % Element-wise multiply: each column of X_fft gets the same H
    % H is [T x 1], X_fft is [T x N] → broadcasting along dim 2
    Y_fft = X_fft .* H;                   % [T x N]

    Y_time = real(ifft(Y_fft, [], 1));    % [T x N] — back to time domain

    % Restore NaN mask: if the original trial was all-NaN, keep it NaN
    all_nan_trials = all(isnan(X), 1);    % [1 x N] logical
    Y = Y_time;
    Y(:, all_nan_trials) = NaN;
end


%% ========================================================================
%%  FUNCTION: plot_wiener_summary
%% ========================================================================
function plot_wiener_summary(times, ...
        ERP_lb_raw, ERP_ps_raw, ...
        ERP_lb_w,   ERP_ps_w, ...
        h, Sdd, Svv, f_axis, Fs, ...
        resid_lb, resid_ps, tag)

    T      = numel(times);
    T_half = floor(T/2) + 1;

    % One-sided spectra (0 .. Nyquist)
    f_one   = f_axis(1:T_half);
    H_gain  = Sdd ./ (Sdd + Svv + eps);
    H_one   = H_gain(1:T_half);
    Sdd_one = Sdd(1:T_half);
    Svv_one = Svv(1:T_half);

    f_plot_max = min(60, Fs/2);
    idx = f_one <= f_plot_max;

    fig = figure('Name', tag, 'Color', 'w', ...
                 'Units', 'normalized', 'Position', [0.05 0.05 0.88 0.82]);

    %% Panel 1 — ERP waveforms
    subplot(2, 2, 1);
    plot(times, ERP_lb_raw, 'Color', [0.6 0.6 0.6],   'LineWidth', 1.2, ...
         'DisplayName', 'Pre (raw)');     hold on;
    plot(times, ERP_ps_raw, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, ...
         'DisplayName', 'Post (raw)');
    plot(times, ERP_lb_w,   'Color', [0.18 0.55 0.34], 'LineWidth', 1.8, ...
         'LineStyle', '--', 'DisplayName', 'Pre (Wiener)');
    plot(times, ERP_ps_w,   'Color', [0.49 0.18 0.56], 'LineWidth', 1.8, ...
         'LineStyle', '--', 'DisplayName', 'Post (Wiener)');
    xline(0, 'k:', 'LineWidth', 0.8);
    yline(0, 'Color', [0.7 0.7 0.7]);
    xlabel('Time (ms)');  ylabel('\muV');
    title('ERP: raw vs Wiener-filtered');
    legend('Location', 'best', 'FontSize', 7);
    grid on;  box off;

    %% Panel 2 — Filter gain H*(f)
    subplot(2, 2, 2);
    fill([f_one(idx); flipud(f_one(idx))], ...
         [H_one(idx); zeros(sum(idx), 1)], ...
         [0.18 0.55 0.34], 'FaceAlpha', 0.15, 'EdgeColor', 'none');  hold on;
    plot(f_one(idx), H_one(idx), 'Color', [0.18 0.55 0.34], 'LineWidth', 1.8);
    yline(0.5, 'k--', 'LineWidth', 0.7);
    text(f_plot_max * 0.6, 0.53, 'H* = 0.5', 'FontSize', 7);
    xlabel('Frequency (Hz)');  ylabel('H*(f)  [0–1]');
    ylim([0 1.05]);
    title('Wiener filter gain  H*(f) = S_{dd}/(S_{dd}+S_{vv})');
    grid on;  box off;

    %% Panel 3 — PSDs + H*(f) overlaid
    subplot(2, 2, 3);
    norm_fac = max(Sdd_one(idx) + Svv_one(idx)) + eps;
    plot(f_one(idx), Sdd_one(idx)/norm_fac, ...
         'Color', [0.49 0.18 0.56], 'LineWidth', 1.5, ...
         'DisplayName', 'S_{dd} (SSEP template)');          hold on;
    plot(f_one(idx), Svv_one(idx)/norm_fac, ...
         'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, ...
         'DisplayName', 'S_{vv} (noise residual)');
    plot(f_one(idx), H_one(idx), ...
         'Color', [0.18 0.55 0.34], 'LineWidth', 1.8, ...
         'LineStyle', '--', 'DisplayName', 'H*(f)');
    xlabel('Frequency (Hz)');  ylabel('Normalised power / gain');
    title('Wiener-Khinchin: PSDs → filter gain');
    legend('Location', 'best', 'FontSize', 7);
    grid on;  box off;

    %% Panel 4 — Impulse response h[n] = ifft(H)
    % fftshift centres the zero-phase response at the middle of the plot.
    subplot(2, 2, 4);
    h_shift = fftshift(h);
    t_h = (-floor(T/2) : floor((T-1)/2)) * (1000/Fs);   % ms, centred
    plot(t_h, h_shift, 'Color', [0.19 0.51 0.74], 'LineWidth', 1.4);
    yline(0, 'Color', [0.7 0.7 0.7]);
    xlabel('Lag (ms)');  ylabel('h[n]');
    title('Impulse response  h = ifft(H*)  [zero-phase]');
    % Mark the main lobe half-width as a rough bandwidth indicator
    [~, pk] = max(abs(h_shift));
    half_amp = h_shift(pk) / 2;
    yline(half_amp, 'r--', 'LineWidth', 0.7);
    grid on;  box off;

    annotation(fig, 'textbox', [0.72 0.01 0.26 0.06], ...
        'String', sprintf('Resid power — Pre: %.4f  |  Post: %.4f  µV²', ...
                          resid_lb, resid_ps), ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', ...
        'FontSize', 7, 'HorizontalAlignment', 'center');

    sgtitle(strrep(tag, '__', ' | '), 'Interpreter', 'none', 'FontSize', 9);
    drawnow;
end