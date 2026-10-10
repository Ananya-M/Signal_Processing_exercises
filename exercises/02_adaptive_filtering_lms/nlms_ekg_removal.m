clc; clear;
%% ========== USER SETTINGS ==============================================

subject_id="Sample_raw_data";
use_nlms = true;   % false = run the same pipeline without NLMS (needed for the QC comparison)

% Sample_data/ lives at the repository root (shared by all exercises).
% Running the whole script: mfilename gives this file's folder -> go up two levels.
% Running a selection: mfilename is empty -> fall back to the current folder.
script_dir = fileparts(mfilename('fullpath'));
data_dir   = fullfile(script_dir, '..', '..', 'Sample_data');
if isempty(script_dir) || ~isfolder(data_dir)
    data_dir = fullfile(pwd, 'Sample_data');
end
if ~isfolder(data_dir)
    error('Sample_data folder not found. Put it in the repository root, or set data_dir by hand.');
end
EEG = pop_loadset('Sample_raw_data.set', data_dir);

bs_start        = -220;
bs_end          = -20;
epoch_interval  = [-0.3 0.7];
Trigger_channels    = {'TRIGGERS4','TRIGGERS8'};
EMG_channels        = {'BIP 01','BIP 02'};
EKG_chan            = {'BIP 04'};
stim_side           = ["L","R"];
selected_EEG_labels={'CREF','Cz','Fz','C3', 'C4','CP5','CP1','CP2','CP6','F4'};
ref='Fz';
EEG.subject_id = subject_id;
channel_labels = {EEG.chanlocs.labels};


%% remove extra chans
EEG = pop_select(EEG, 'nochannel', {'X-AXIS','Y-AXIS','Z-AXIS','TRIGGERS1','TRIGGERS2', ...
    'TRIGGERS3','TRIGGERS5','TRIGGERS6','TRIGGERS7','TRIGGERS9',...
    'TRIGGERS10','TRIGGERS11','TRIGGERS12','TRIGGERS13','TRIGGERS14','TRIGGERS15','TRIGGERS16', ...
    'STATUS','COUNTER','Pz','F3','FC1','P3','FC2','P4'});
%%
selected_EEG_labels = selected_EEG_labels(ismember(selected_EEG_labels, channel_labels));
selected_EEG_channels = find(ismember({EEG.chanlocs.labels}, selected_EEG_labels));

%% import trigger channel
trig_chan = find(ismember({EEG.chanlocs.labels}, Trigger_channels));
all_events = EEG.event; % Store existing events
for i = 1:length(trig_chan)
    EEG_tmp = pop_chanevent(EEG, trig_chan(i), 'edge', 'leading', 'delchan', 'off'); % Extract events

    if ~isempty(EEG_tmp.event) % Ensure new events exist
        new_events = EEG_tmp.event;

        % Assign unique event types (optional)
        for j = 1:length(new_events)
            new_events(j).type = ['Trigger' num2str(trig_chan(i))]; % Label events per channel
        end

        all_events = [all_events, new_events];  % Append new events to existing ones
    end
end
EEG.event = all_events; % Save all events back to EEG
EEG = eeg_checkset(EEG); % Update EEGLAB structure

%%
trigs=unique({(EEG.event.type);});
nums = cellfun(@(s) str2double(regexp(s,'\d+','match','once')), trigs);
[~, idx] = min(nums);
left_trig = trigs{idx};

%% high pass for cleanline
[EEG, com, b] = pop_eegfiltnew(EEG, 'locutoff', 1,'channels', selected_EEG_channels);
EEG = eeg_checkset(EEG);

%% cleanline
EEG = eeg_checkset(EEG);
figure;
pop_spectopo(EEG, 1, [], 'EEG', 'freqrange', [0.1 100]);
EEG = pop_cleanline(EEG, 'Bandwidth', 2, ... % Taper bandwidth = 2 Hz
    'ChanCompIndices', selected_EEG_channels, ...
    'SignalType', 'Channels', ...
    'ComputeSpectralPower', true, ...
    'LineFrequencies', [60], ... % Removes 60 Hz line noise
    'NormalizeSpectrum', false, ...
    'LineAlpha', 0.01, ... % Significance threshold for detecting line noise
    'PaddingFactor', 2, ...
    'PlotFigures', false, ...
    'ScanForLines', false, ...
    'SmoothingFactor', 100, ...
    'VerbosityLevel', 1, ...
    'SlidingWinLength', 4, ...
    'SlidingWinStep', 2);
figure;
pop_spectopo(EEG, 1, [], 'EEG', 'freqrange', [1 100]);


%% bandpass
EEG = eeg_checkset(EEG);
selected_BIP_channels=find(ismember({EEG.chanlocs.labels}, EMG_channels));
[EEG, com, b] = pop_eegfiltnew(EEG, 'locutoff', 1,'hicutoff',50,'channels', selected_EEG_channels);
[EEG, com, b] = pop_eegfiltnew(EEG, 'locutoff', 1,'channels', selected_BIP_channels);


%% remove EKG artifact (NLMS, continuous data, before epoching)
stim_types = unique({EEG.event.type});
stim_types = stim_types(startsWith(stim_types, 'Trigger'));            % both trigger channels
stim_samples = [EEG.event(ismember({EEG.event.type}, stim_types)).latency];   % continuous, samples
if use_nlms
    [EEG, EKG_NLMS] = remove_ekg_nlms(EEG, EKG_chan{1}, selected_EEG_labels, ...
        'FilterDurationMs', 100, 'TapSpacingMs', 2, 'Mu', 0.01, ...
        'BlankLatencies', stim_samples, 'BlankWindowMs', [-10 15]);
else
    % no filtering: keep only the mean-removed EKG reference for the phase-locking check
    ekg_idx = find(strcmp({EEG.chanlocs.labels}, EKG_chan{1}));
    ref_raw = double(EEG.data(ekg_idx,:));
    ref_raw = ref_raw - mean(ref_raw);
    EKG_NLMS = struct('ekg_ref_used', ref_raw, 'ekg_ref_raw', ref_raw);
end
EEG.EKG_NLMS = EKG_NLMS;

%% plotting: before / after NLMS (first 240 s)
if use_nlms
    plot_chan = 'CP2';
    ch_idx = find(strcmp({EEG.chanlocs.labels}, plot_chan));   % index into EEG.data, not into selected_EEG_labels
    n_plot = min(240*EEG.srate, EEG.pnts);
    figure;
    plot(EEG.times(1:n_plot), EEG.data_before_ekg_removal(ch_idx,1:n_plot), 'DisplayName','before'); hold on;
    plot(EEG.times(1:n_plot), EEG.data(ch_idx,1:n_plot), 'DisplayName','after');
    xlabel('Time (ms)'); ylabel('EEG');
    legend;
    title(sprintf('EEG at %s before and after NLMS filtering', plot_chan));
    plot_field = matlab.lang.makeValidName(plot_chan);
    figure;
    plot(EKG_NLMS.(plot_field).block_mse, '-o');
    xlabel('Time (s, 1s blocks)'); ylabel('Mean error power');
    title(sprintf('EKG-NLMS convergence — %s', plot_chan));
end

%% ===== CAPTURE CONTINUOUS STIM TIMES + PHASE-LOCKING CHECK =====
% MUST run BEFORE pop_epoch. After epoching, EEG.event.latency is re-indexed
% into the concatenated epoched data (roughly (k-1)*pnts + time-zero offset),
% which no longer lines up with the continuous R-peak times.
fs = EEG.srate;

trig_types   = unique({EEG.event.type});
right_trigs  = setdiff(trig_types(startsWith(trig_types,'Trigger')), left_trig);
stim_left_s  = ([EEG.event(strcmp({EEG.event.type}, left_trig)).latency]    - 1) / fs;
stim_right_s = ([EEG.event(ismember({EEG.event.type}, right_trigs)).latency] - 1) / fs;

% keep them for later use anywhere in the pipeline
EEG.stim_latency_cont_s = struct('left', stim_left_s, 'right', stim_right_s);

% sanity checks: continuous timeline looks like real stim spacing, not epoch spacing
rec_dur_s = (EEG.pnts - 1) / fs;
fprintf('Recording duration: %.1f s\n', rec_dur_s);
fprintf('Left  stims: n=%d, median inter-stim interval = %.3f s, last stim at %.1f s\n', ...
    numel(stim_left_s),  median(diff(stim_left_s)),  stim_left_s(end));
fprintf('Right stims: n=%d, median inter-stim interval = %.3f s, last stim at %.1f s\n', ...
    numel(stim_right_s), median(diff(stim_right_s)), stim_right_s(end));
% If the median interval is a clean multiple of epoch length (e.g. 1.000 s) the
% latencies are epoched indices and this check is invalid.

% R-peak detection on a QRS-band copy (does NOT touch the NLMS reference or data)
ekg = double(EEG.EKG_NLMS.ekg_ref_raw);   % unblanked, so R-peaks inside blanked windows are not lost
[bQ,aQ] = butter(2, [5 25]/(fs/2), 'bandpass');
ekg_d = filtfilt(bQ, aQ, ekg);
if abs(min(ekg_d)) > abs(max(ekg_d)), ekg_d = -ekg_d; end   % force R-peaks positive
[~, r_locs] = findpeaks(ekg_d, 'MinPeakHeight', 3*std(ekg_d), ...
    'MinPeakDistance', round(0.4*fs));
r_peak_times = (r_locs - 1) / fs;

% visual check of detection on first 20 s
figure;
seg = 1:min(20*fs, numel(ekg_d));
plot(seg/fs, ekg_d(seg)); hold on;
rp = r_locs(r_locs <= seg(end));
plot(rp/fs, ekg_d(rp), 'rv', 'MarkerFaceColor','r');
xlabel('Time (s)'); title('R-peak detection check (first 20 s)');

res_left  = ekg_stim_phaselock(r_peak_times, stim_left_s,  'Left trigger');
res_right = ekg_stim_phaselock(r_peak_times, stim_right_s, 'Right trigger');
EEG.phaselock = struct('left', res_left, 'right', res_right);

%% Phase vs trial number (looks for beat patterns between stim rate and heart rate)
if ~exist('res_left','var') || ~exist('res_right','var')
    res_left  = EEG.phaselock.left;
    res_right = EEG.phaselock.right;
end
figure; plot(res_right.phase, '.-'); xlabel('Trial #'); ylabel('Cardiac phase');
hold on; plot(res_left.phase, '.-');
legend('Right trigger','Left trigger');


%% epoch
EEG = eeg_checkset(EEG);
EEG = pop_epoch(EEG, stim_types, epoch_interval);   % stim triggers only, not every event type


%% remove baseline
EEG.raw_data = EEG.data;
EEG = eeg_checkset(EEG);
[EEG,V] = remove_epoch_baseline(EEG, bs_start, bs_end);

plot_early_ERP(EEG,'raw_data');
figure;
plot_early_ERP(EEG,'data');


%% change ref

EEG = eeg_checkset(EEG);
ref_chan = find(ismember({EEG.chanlocs.labels}, ref));
EEG = pop_reref(EEG, ref_chan);
plot_early_ERP(EEG,'data');

%% remove bad trials

EEG = eeg_checkset(EEG);

channel_labels = {EEG.chanlocs.labels};
ignore_chan = find(~cellfun('isempty', ...
    regexp(channel_labels,'(trigger|bip)','ignorecase')));

% ---- target rejection range ----
minPct = 0.05;
maxPct = 0.10;
nTrials = EEG.trials;
rmepochs = [];
startprob_list = [ 3 4 4.5 ];

pct_list   = NaN(size(startprob_list));
epochs_cell = cell(size(startprob_list));

for k = 1:numel(startprob_list)
    sp = startprob_list(k);
    EEG_tmp = EEG;

    [~, rmepochs] = artifact_detection(EEG_tmp, sp, ignore_chan);

    pctRejected = numel(rmepochs) / nTrials;

    fprintf('startprob %.1f → %.2f%% rejected\n', sp, 100*pctRejected);

    pct_list(k)    = pctRejected;
    epochs_cell{k} = rmepochs;
end

valid_idx = find(pct_list >= minPct & pct_list <= maxPct);

if isempty(valid_idx)
    [~, best_local] = min(abs(pct_list - maxPct));
    best_idx = best_local;
    warning('No startprob produced rejection between %.0f–%.0f%%', ...
        minPct*100, maxPct*100);
else
    [~, best_local] = min(abs(pct_list(valid_idx) - maxPct));
    best_idx = valid_idx(best_local);
end

startprob_used = startprob_list(best_idx);
rmepochs       = epochs_cell{best_idx};
pctRej         = pct_list(best_idx);

fprintf('\n✅ Selected startprob %.1f → %.2f%% rejected\n', ...
    startprob_used, 100*pctRej);

EEG.startprob     = startprob_used;
EEG.badtrial_pct  = pctRej;

%%
EEG = eeg_checkset(EEG);
EEG.data_without_badtrials = EEG.data;
EEG.badtrials=rmepochs;

%%
right_chans={'C4','CP2','CP6'};
left_chans={'C3','CP1','CP5'};
%% ===================== BUILD HDS DATA =====================

EEG.hds_data = struct();
right_output_names = {'rightmn_data_without_badtrials'};
left_output_names  = {'leftmn_data_without_badtrials'};
data_fields        = {'data_without_badtrials'};

combField_right = ['combined_' strjoin(right_chans,'')];
combField_left  = ['combined_' strjoin(left_chans,'')];

%% LEFT Trigger
for d = 1:numel(data_fields)
    current_field = data_fields{d};
    [left_cond,~] = splitTrials(EEG, left_trig);
    EEG.hds_data.(left_output_names{d}) = struct();

    for k = 1:length(selected_EEG_labels)
        chanLabel = selected_EEG_labels{k};
        chanField = matlab.lang.makeValidName(chanLabel);
        chan_idx = find(strcmp({EEG.chanlocs.labels}, chanLabel));
        tmp = EEG.(current_field)(chan_idx,:,left_cond);
        tmp = permute(tmp,[2 3 1]);   % time x trials
        EEG.hds_data.(left_output_names{d}).(chanField) = tmp;
    end

    EEG.hds_data.(left_output_names{d}) = ...
        addCombinedChannels(EEG.hds_data.(left_output_names{d}), right_chans, combField_right);
end

%% RIGHT Trigger
for d = 1:numel(data_fields)
    current_field = data_fields{d};
    [~,right_cond] = splitTrials(EEG,left_trig);
    EEG.hds_data.(right_output_names{d}) = struct();

    for k = 1:length(selected_EEG_labels)
        chanLabel = selected_EEG_labels{k};
        chanField = matlab.lang.makeValidName(chanLabel);
        chan_idx = find(strcmp({EEG.chanlocs.labels}, chanLabel));
        tmp = EEG.(current_field)(chan_idx,:,right_cond);
        tmp = permute(tmp,[2 3 1]);
        EEG.hds_data.(right_output_names{d}).(chanField) = tmp;
    end

    EEG.hds_data.(right_output_names{d}) = ...
        addCombinedChannels(EEG.hds_data.(right_output_names{d}), left_chans, combField_left);
end

%% ====================== SAVE & VISUALIZE ======================

folder_save = fullfile(data_dir, 'processed_data_LMS');
if ~exist(folder_save, 'dir')
    mkdir(folder_save);
end
figure; plot_early_ERP(EEG, 'data_without_badtrials');
if use_nlms
    save_filename = subject_id + "_with_LMS.set";
else
    save_filename = subject_id + "_without_LMS.set";
end
EEG = pop_saveset(EEG, 'filename', char(save_filename), 'filepath', folder_save);
disp('✅ Processing and save completed successfully.');


%% QC plotting
% Needs both runs: run the script once with use_nlms = true and once with
% use_nlms = false. Each run saves its own file in folder_save.
file_without = char(subject_id + "_without_LMS.set");
file_with    = char(subject_id + "_with_LMS.set");
if ~isfile(fullfile(folder_save, file_without)) || ~isfile(fullfile(folder_save, file_with))
    fprintf('QC skipped: run the script with use_nlms = true and use_nlms = false first.\n');
    return
end
EEG1 = pop_loadset(file_without, folder_save);   % no NLMS
EEG2 = pop_loadset(file_with,    folder_save);   % NLMS
ERP_1 = squeeze(nanmean(EEG1.hds_data.leftmn_data_without_badtrials.CP6(:, :), 2));
ERP_2 = squeeze(nanmean(EEG2.hds_data.leftmn_data_without_badtrials.CP6(:, :), 2));

figure;
plot(EEG.times, ERP_1, 'LineWidth', 2, 'DisplayName','ERP without LMS filtering'); hold on;
plot(EEG.times, ERP_2, 'LineWidth', 2, 'DisplayName','ERP with LMS filtering');
xline(0, '--b', 'LineWidth', 2, 'DisplayName','Stimulus'); hold off
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title(sprintf('ERP comparison of LMS filtering on (%s)', 'CP6'));
xlim([-300 300]);
legend('show', 'Location', 'bestoutside');
grid on;


%% some additional checks - check where exactly are the differences

A = EEG1.hds_data.leftmn_data_without_badtrials.CP6;   % time x trials, no NLMS
B = EEG2.hds_data.leftmn_data_without_badtrials.CP6;   % time x trials, with NLMS
fprintf('trials: %d (no NLMS) vs %d (NLMS)\n', size(A,2), size(B,2));
nT = min(size(A,2), size(B,2));
if size(A,2) ~= size(B,2)
    warning('Trial counts differ; pairing the first %d trials. Rerun both versions with the same script so the trial sets match.', nT);
end
D = B(:,1:nT) - A(:,1:nT);                           % what NLMS changed, per trial
t = EEG2.times;  m = mean(D,2);  sem = std(D,0,2)/sqrt(size(D,2));
figure; plot(t,m); hold on; plot(t,m+2*sem,':k', t,m-2*sem,':k');
xline(0,'--'); xlim([-300 300]); xlabel('Time (ms)'); ylabel('NLMS minus no-NLMS (\muV)');

%% how does different u affect
ch  = find(strcmp({EEG2.chanlocs.labels}, 'CP6'));
d   = double(EEG2.data_before_ekg_removal(ch,:));
ref_sig = double(EEG2.EKG_NLMS.ekg_ref_raw);      % unblanked, mean-removed (named so it doesn't overwrite ref = 'Fz')
lags = EEG2.EKG_NLMS.lags;  fs = EEG.srate;

e_real  = nlms_1ch(d, ref_sig, lags, 0.01);                          % pipeline settings
e_shuf  = nlms_1ch(d, circshift(ref_sig, round(60*fs)), lags, 0.01); % same reference, wrong timing
e_lowmu = nlms_1ch(d, ref_sig, lags, 0.001);                         % 10x smaller step

stim = round(EEG2.stim_latency_cont_s.left*fs) + 1;
pre = round(0.3*fs);  post = round(0.3*fs);
stim = stim(stim-pre >= 1 & stim+post <= numel(d));
t = (-pre:post)/fs*1000;  bl = t >= -220 & t <= -20;
erp = @(x) mean(cell2mat(arrayfun(@(i) x(i-pre:i+post) - mean(x(i-pre+find(bl)-1)), ...
              stim(:), 'UniformOutput', false)), 1);

figure; hold on;
plot(t, erp(d),       'k', 'LineWidth',2, 'DisplayName','no NLMS');
plot(t, erp(e_real),  'DisplayName','NLMS, real ref, \mu=0.01');
plot(t, erp(e_shuf),  'DisplayName','NLMS, shuffled ref');
plot(t, erp(e_lowmu), 'DisplayName','NLMS, real ref, \mu=0.001');
xline(0,'--'); xlim([-100 200]); legend; xlabel('Time (ms)'); ylabel('\muV');

%% null test: NLMS effect at random pseudo-stim times vs. real stims
% reuses d, e_real, fs, stim from the previous section
pre = round(0.3*fs); post = round(0.3*fs);
t = (-pre:post)/fs*1000; bl = t >= -220 & t <= -20; win = t >= 0 & t <= 100;
n = numel(stim); K = 300;
valid = (pre+1):(numel(d)-post);

rms_d = zeros(K,1); rms_e = zeros(K,1); rms_diff = zeros(K,1);
for k = 1:K
    ps = valid(randperm(numel(valid), n));       % random pseudo-stims, same count
    A = erp_at(d,      ps, pre, post, bl);
    B = erp_at(e_real, ps, pre, post, bl);
    rms_d(k) = sqrt(mean(A(win).^2));
    rms_e(k) = sqrt(mean(B(win).^2));
    rms_diff(k) = sqrt(mean((B(win)-A(win)).^2));
end
A = erp_at(d, stim, pre, post, bl);  B = erp_at(e_real, stim, pre, post, bl);
real_diff = sqrt(mean((B(win)-A(win)).^2));

fprintf('Null ERP RMS (0-100 ms): no NLMS %.2f, NLMS %.2f uV\n', mean(rms_d), mean(rms_e));
fprintf('NLMS-minus-raw RMS: null median %.2f (95th pct %.2f) vs real stims %.2f uV\n', ...
    median(rms_diff), prctile(rms_diff,95), real_diff);


%% ====================== HELPER FUNCTIONS ======================

function e = nlms_1ch(d, ref, lags, mu)
% Single-channel NLMS, same update as remove_ekg_nlms (used by the QC control runs).
    M = lags(end) + 1;  N = numel(d);
    ref_padded = [zeros(1,M) ref];
    w = zeros(numel(lags),1);  e = zeros(1,N);
    for n = 1:N
        x = ref_padded(M + n - lags)';
        e(n) = d(n) - w'*x;
        w = w + (mu/(x'*x + 1e-6)) * e(n) * x;
    end
end


function m = erp_at(x, idx, pre, post, bl)
% Baseline-corrected average of x around each index in idx (used by the null test).
    S = cell2mat(arrayfun(@(i) x(i-pre:i+post), idx(:), 'UniformOutput', false));
    S = S - mean(S(:,bl), 2);
    m = mean(S, 1);
end


function res = ekg_stim_phaselock(r_peak_times, stim_times, label)
% Cardiac phase at each stim onset and a Rayleigh test for non-uniformity.
% Phase = (stim - previous R) / (next R - previous R), in [0,1). Using the
% local RR interval (not a global median) avoids smearing from heart-rate
% variability, and one-sided "time since last R" avoids the RR/2 ceiling of a
% nearest-neighbour |distance|, which fakes clustering under the null.
    median_RR = median(diff(r_peak_times));

    phase = nan(size(stim_times));
    for i = 1:numel(stim_times)
        prev_r = find(r_peak_times <  stim_times(i), 1, 'last');
        next_r = find(r_peak_times >= stim_times(i), 1, 'first');
        if isempty(prev_r) || isempty(next_r), continue; end
        rr_local = r_peak_times(next_r) - r_peak_times(prev_r);
        if rr_local > 1.8*median_RR, continue; end          % missed beat -> phase unreliable
        phase(i) = (stim_times(i) - r_peak_times(prev_r)) / rr_local;
    end
    phase = phase(~isnan(phase));

    theta = 2*pi*phase;
    n  = numel(theta);
    C  = mean(exp(1i*theta));
    R  = abs(C);
    z  = n * R^2;
    Rn = n * R;
    % Zar (1999) Rayleigh p-value approximation, valid for small and large z
    p  = exp(sqrt(1 + 4*n + 4*(n^2 - Rn^2)) - (1 + 2*n));
    p  = min(max(p, 0), 1);

    fprintf('[%s] median RR = %.3f s | n = %d | R = %.3f | z = %.3f | p = %.4g | mean phase = %.2f cycles\n', ...
        label, median_RR, n, R, z, p, mod(angle(C)/(2*pi), 1));

    figure;
    subplot(1,2,1); histogram(phase, 10, 'BinLimits',[0 1]);
    xlabel('Cardiac phase at stim (0 = R-peak, 1 = next R-peak)'); ylabel('Count');
    title(sprintf('%s: should be flat if independent', label));
    subplot(1,2,2); polarhistogram(theta, 12);
    title(sprintf('R = %.3f, p = %.3g, n = %d', R, p, n));

    res = struct('R',R, 'z',z, 'p',p, 'n',n, 'median_RR',median_RR, 'phase',phase);
end


function dataStruct = addCombinedChannels(dataStruct, chanSet, combField)
    combined = combineChannels(dataStruct, chanSet);
    dataStruct.(combField) = combined;
end


function [eeg, rmepochs] = artifact_detection(eeg, x, ignore_chan)
[~, rmepochs] = pop_autorej(eeg, ...
    'nogui','on', ...
    'electrodes', setdiff(1:eeg.nbchan, ignore_chan), ...
    'threshold', 100, ...
    'startprob', x);
end


% ---------- Split trials ----------
function [left_trials,right_trials] = splitTrials(EEG,left_trig)
event_types  = {EEG.event.type};
event_trials = [EEG.event.epoch];
left_trials  = event_trials(strcmp(event_types,left_trig));
right_trials = setdiff(event_trials,left_trials);
end


function combined = combineChannels(chanData, chanSet)
    allData = cell(1,numel(chanSet));
    for ch = 1:numel(chanSet)
        dat = chanData.(chanSet{ch});
        if size(dat,1)==1
            dat = dat';
        end
        allData{ch} = dat;
    end
    combined = cat(2, allData{:});
end


function [EEG, EKG_NLMS] = remove_ekg_nlms(EEG, ekg_label, target_labels, varargin)
% REMOVE_EKG_NLMS  Continuous-data NLMS cancellation of the ELECTRICAL
% QRS artifact from EEG channels (BCG/mechanical component NOT targeted).
% Run on continuous data, BEFORE epoching.
%
% Name-value options:
%   'FilterDurationMs' - filter memory in ms (default 100)
%   'TapSpacingMs'     - currently ignored: tap spacing fixed at 2 SAMPLES
%   'Mu'               - NLMS step size (default 0.01)
%   'Epsilon'          - regularizer for near-silent reference (default 1e-6)
%   'RefBandpassHz'    - bandpass on the EKG reference ONLY (default [] = skip)

p = inputParser;
addParameter(p,'FilterDurationMs',100);
addParameter(p,'TapSpacingMs',2);
addParameter(p,'Mu',0.01);
addParameter(p,'Epsilon',1e-6);
addParameter(p,'RefBandpassHz',[]);
addParameter(p,'BlankLatencies',[]);        % continuous stim latencies, in samples
addParameter(p,'BlankWindowMs',[-5 15]);    % window around each stim to blank, ms
parse(p,varargin{:});
opt = p.Results;

fs = EEG.srate;
tap_spacing = 2;
n_taps      = round(opt.FilterDurationMs/1000*fs / tap_spacing);
lags        = (0:n_taps-1) * tap_spacing;
M           = lags(end) + 1;

ekg_idx = find(strcmp({EEG.chanlocs.labels}, ekg_label));
if isempty(ekg_idx)
    error('EKG channel "%s" not found in EEG.chanlocs.', ekg_label);
end
target_idx = find(ismember({EEG.chanlocs.labels}, target_labels));
if isempty(target_idx)
    error('None of the requested target_labels were found.');
end

ref = double(EEG.data(ekg_idx,:));
ref = ref - mean(ref);
if ~isempty(opt.RefBandpassHz)
    [bB,bA] = butter(4, opt.RefBandpassHz/(fs/2), 'bandpass');
    ref = filtfilt(bB, bA, ref);
end

ref_raw = ref;   % unblanked copy, used for R-peak detection

% Blank the reference around every stim (linear interpolation across the gap) so it
% carries no trigger-locked marker the filter could learn and subtract from the EEG.
if ~isempty(opt.BlankLatencies)
    a = round(opt.BlankWindowMs(1)/1000*fs);
    b = round(opt.BlankWindowMs(2)/1000*fs);
    nref = numel(ref);
    for k = 1:numel(opt.BlankLatencies)
        i1 = round(opt.BlankLatencies(k)) + a;
        i2 = round(opt.BlankLatencies(k)) + b;
        if i1 < 2 || i2 > nref-1, continue; end
        ref(i1:i2) = interp1([i1-1 i2+1], [ref(i1-1) ref(i2+1)], i1:i2);
    end
end

N = length(ref);
ref_padded = [zeros(1, M), ref];

EEG.data_before_ekg_removal = EEG.data;

EKG_NLMS = struct('params', opt, 'lags', lags, 'fs', fs, ...
    'ekg_ref_used', ref, 'ekg_ref_raw', ref_raw);   % used = blanked (fed to NLMS), raw = unblanked

for c = 1:numel(target_idx)
    chan_i   = target_idx(c);
    chanName = matlab.lang.makeValidName(EEG.chanlocs(chan_i).labels);
    d = double(EEG.data(chan_i,:));

    w = zeros(n_taps,1);
    e = zeros(1,N);
    err_energy = zeros(1,N);

    for n = 1:N
        x_n = ref_padded(M + n - lags)';
        y_n = w' * x_n;
        e(n) = d(n) - y_n;
        norm_factor = (x_n' * x_n) + opt.Epsilon;
        w = w + (opt.Mu / norm_factor) * e(n) * x_n;
        err_energy(n) = e(n)^2;
    end

    EEG.data(chan_i,:) = e;

    EKG_NLMS.(chanName).final_weights = w;
    block = fs;
    nblocks = floor(N/block);
    EKG_NLMS.(chanName).block_mse = arrayfun(@(b) ...
        mean(err_energy((b-1)*block+1:b*block)), 1:nblocks);
end

end
